from __future__ import annotations

import logging
import os
import re
from pathlib import Path
from typing import Dict

import anndata as ad
import pandas as pd
from sqlalchemy.orm import Session

from histocat.core.acquisition import service as acquisition_service
from histocat.core.dataset import service as dataset_service
from histocat.core.dataset.dto import DatasetCreateDto, DatasetUpdateDto
from histocat.core.dataset.models import CELL_FILENAME, DatasetModel
from histocat.core.notifier import Message
from histocat.core.project import service as project_service
from histocat.core.redis_manager import UPDATES_CHANNEL_NAME, redis_manager
from histocat.core.slide import service as slide_service
from histocat.worker.io.utils import CELL_CSV_FILENAME, IMAGE_CSV_FILENAME, copy_file

logger = logging.getLogger(__name__)


PANEL_CSV_FILE = "panel.csv"


def import_dataset(db: Session, root_folder: Path, cell_csv_filename: str, project_id: int):
    """Import dataset from the folder compatible with 'cpout' IMC pipeline folders."""

    project = project_service.get(db, id=project_id)
    if not project:
        logger.warning(f"Cannot import dataset: project [id: {project_id}] does not exist.")
        return

    create_params = DatasetCreateDto(project_id=project_id, origin="ImcSegmentationPipelineV2", status="pending")
    dataset = dataset_service.create(db, params=create_params)
    meta = {}

    src_folder = Path(cell_csv_filename).parent
    dst_folder = Path(dataset.location)
    os.makedirs(dst_folder, exist_ok=True)

    image_df = _import_image_csv(src_folder)

    masks = {}
    image_number_to_acquisition_id = {}
    image_number_to_scaling = {}
    for index, row in image_df.iterrows():
        mask_meta = _import_mask(db, src_folder, row, dataset)
        if mask_meta is not None:
            acquisition_id = mask_meta.get("acquisition").get("id")
            masks[acquisition_id] = mask_meta
            image_number = mask_meta.get("image_number")
            image_number_to_acquisition_id[image_number] = acquisition_id

            scaling = int(row["Scaling_FullStack"])
            image_number_to_scaling[image_number] = scaling
    meta["masks"] = masks

    # Import panel data: { Metal Tag : channel number }
    channel_order, channel_mass = _import_panel(os.path.join(src_folder, PANEL_CSV_FILE))

    # Convert cell.csv to AnnData file format
    cell_df = _import_cell_csv(
        src_folder, dst_folder, image_number_to_acquisition_id, image_number_to_scaling, channel_order
    )

    # Register neighbors columns
    neighbors_cols = [col.split("_")[1] for col in cell_df.columns if "Neighbors_" in col]
    meta["columns"] = {"neighbors": neighbors_cols}

    acquisition_ids = sorted(list(masks.keys()))
    channels_sorted = sorted(channel_mass.items(), key=lambda kv: kv[1])
    channels = [c[0] for c in channels_sorted]

    update_params = DatasetUpdateDto(
        name=f"Dataset {dataset.id}", status="ready", acquisition_ids=acquisition_ids, channels=channels, meta=meta
    )
    dataset = dataset_service.update(db, item=dataset, params=update_params)
    redis_manager.publish(UPDATES_CHANNEL_NAME, Message(project_id, "dataset_imported"))


def _import_image_csv(src_folder: Path):
    src_uri = src_folder / IMAGE_CSV_FILENAME

    if not src_uri.exists():
        return None

    df = pd.read_csv(src_uri)
    return df


def _import_panel(src: str):
    df = pd.read_csv(src)
    # Map Metal Tag to its order number and mass
    channel_order = {}
    channel_mass = {}
    for _, row in df.iterrows():
        channel_order[row.loc["Metal Tag"]] = int(row.loc["channel"])
        channel_mass[row.loc["Metal Tag"]] = int(row.loc["mass"])
    return channel_order, channel_mass


def _import_cell_csv(
    src_folder: Path,
    dst_folder: Path,
    image_number_to_acquisition_id: Dict[int, int],
    image_number_to_scaling: Dict[int, int],
    channel_order: Dict[str, int],
):
    src_uri = src_folder / CELL_CSV_FILENAME

    if not src_uri.exists():
        return None

    df = pd.read_csv(src_uri)
    df.index.astype(str, copy=False)

    obs = pd.DataFrame(index=df.index)
    obs["CellId"] = df.index
    obs["AcquisitionId"] = df["ImageNumber"]
    obs["AcquisitionId"].replace(image_number_to_acquisition_id, inplace=True)
    obs["ImageNumber"] = df["ImageNumber"]
    obs["ObjectNumber"] = df["ObjectNumber"]
    obs["CentroidX"] = df["Location_Center_X"]
    obs["CentroidY"] = df["Location_Center_Y"]

    neighbors_cols = [col for col in df.columns if "Neighbors_" in col]
    for col in neighbors_cols:
        col_name = col.split("_")[1]
        if col_name:
            obs[col_name] = df[col]
            if col_name == "NumberOfNeighbors":
                obs[col_name] = obs[col_name].astype("int64")

    var_names = []
    x_df = pd.DataFrame()
    for key, value in channel_order.items():
        # TODO: check intensity multiplier
        x_df[key] = df[f"Intensity_MeanIntensity_FullStackFiltered_c{value}"] * 2 ** 16
        var_names.append(key)
    var = pd.DataFrame(index=var_names)
    var["Channel"] = var.index
    X_counts = x_df.to_numpy()

    adata = ad.AnnData(X_counts, obs=obs, var=var, dtype="float32")
    dst_uri = dst_folder / CELL_FILENAME
    adata.write_h5ad(dst_uri)

    return df


def _import_mask(db: Session, src_folder: Path, row: pd.Series, dataset: DatasetModel):
    filename = row["FileName_CellImage"]
    image_number = row["ImageNumber"]
    uri = src_folder / "masks" / filename

    p = re.compile("(?P<Name>.*)_s(?P<SlideID>[0-9]+)_a(?P<AcquisitionID>[0-9]+)_ac.*")
    slide_name, slide_origin_id, acquisition_origin_id = p.findall(filename)[0]

    slide = slide_service.get_by_name(db, project_id=dataset.project_id, name=slide_name)
    if slide is None:
        return None
    acquisition = acquisition_service.get_by_origin_id(db, slide_id=slide.id, origin_id=acquisition_origin_id)

    location = copy_file(uri, dataset.location)

    meta = {
        "image_number": image_number,
        "location": location,
        "slide": {"id": slide.id, "origin_id": slide.origin_id},
        "acquisition": {"id": acquisition.id, "origin_id": acquisition.origin_id},
    }
    return meta
