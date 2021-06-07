import logging
import os
import re
from pathlib import Path
from typing import Dict, Union, List

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


def import_dataset(db: Session, input_folder: Path, project_id: int):
    """Import dataset from the folder compatible with 'steinbock' format."""

    project = project_service.get(db, id=project_id)
    if not project:
        logger.warning(f"Cannot import dataset: project [id: {project_id}] does not exist.")
        return

    src_folder = None
    for path in input_folder.rglob(PANEL_CSV_FILE):
        src_folder = path.parent
        break

    if src_folder is None:
        return

    create_params = DatasetCreateDto(project_id=project_id, origin="steinbock", status="pending")
    dataset = dataset_service.create(db, params=create_params)

    # Metadata dictionary
    meta = {}

    dst_folder = Path(dataset.location)
    os.makedirs(dst_folder, exist_ok=True)

    # Import panel data: { Metal Tag : channel number }
    channel_order = _import_panel(os.path.join(src_folder, PANEL_CSV_FILE))

    masks = {}
    acquisition_id_mapping = {}
    for mask_file in sorted(Path(src_folder / "cell_masks").rglob("*.tiff")):
        result = _import_mask(db, mask_file, dataset)
        if result is not None:
            mask_meta, slide_name, acquisition_origin_id = result
            acquisition_id = mask_meta.get("acquisition").get("id")
            image_number = mask_meta.get("acquisition").get("origin_id")
            masks[acquisition_id] = mask_meta
            acquisition_id_mapping[f"{slide_name}_{image_number}"] = acquisition_id
    meta["masks"] = masks

    regionprops_df = pd.DataFrame()
    for regionprops_file in sorted(Path(src_folder / "cell_regionprops").rglob("*.csv")):
        slide_name, acquisition_origin_id, df = _import_regionprops(regionprops_file, acquisition_id_mapping)
        regionprops_df = regionprops_df.append(df)
    regionprops_df.reset_index(inplace=True, drop=True)
    regionprops_df["CellId"] = regionprops_df.index

    obs = pd.DataFrame(index=regionprops_df.index)
    obs["CellId"] = regionprops_df.index
    obs["AcquisitionId"] = regionprops_df["AcquisitionId"]
    obs["ImageNumber"] = regionprops_df["ImageNumber"]
    obs["ObjectNumber"] = regionprops_df["ObjectNumber"]
    obs["CentroidX"] = regionprops_df["CentroidX"]
    obs["CentroidY"] = regionprops_df["CentroidY"]

    print(obs)

    var_names = []
    x_df = pd.DataFrame()
    for key, value in channel_order.items():
        # TODO: check intensity multiplier
        x_df[key] = range(348) # df[f"Intensity_MeanIntensity_FullStackFiltered_c{value}"] * 2 ** 16
        var_names.append(key)
    var = pd.DataFrame(index=var_names)
    var["Channel"] = var.index
    X_counts = x_df.to_numpy()

    adata = ad.AnnData(X_counts, obs=obs, var=var, dtype="float32")
    dst_uri = dst_folder / CELL_FILENAME
    adata.write_h5ad(dst_uri)

    # # Convert cell.csv to AnnData file format
    # cell_df = _import_cell_csv(src_folder, dst_folder, image_number_to_acquisition_id, channel_order)

    acquisition_ids = sorted(list(masks.keys()))
    channels = [c[0] for c in channel_order]

    update_params = DatasetUpdateDto(
        name=f"Dataset {dataset.id}", status="ready", acquisition_ids=acquisition_ids, channels=channels, meta=meta
    )
    dataset = dataset_service.update(db, item=dataset, params=update_params)
    redis_manager.publish(UPDATES_CHANNEL_NAME, Message(project_id, "dataset_imported"))


def _import_panel(path: str):
    panel_df = pd.read_csv(path)
    # Map Metal Tag to its order number
    channel_order = dict(
        [(metal_name, int(index)) for metal_name, index in zip(panel_df.channel, panel_df.index)]
    )
    return channel_order


def _import_mask(db: Session, filepath: Path, dataset: DatasetModel):
    p = re.compile("(?P<Name>.*)_(?P<AcquisitionID>[0-9]+).tiff")
    slide_name, acquisition_origin_id = p.findall(filepath.name)[0]

    slide = slide_service.get_by_name(db, project_id=dataset.project_id, name=slide_name)
    if slide is None:
        return None
    acquisition = acquisition_service.get_by_origin_id(db, slide_id=slide.id, origin_id=acquisition_origin_id)

    location = copy_file(str(filepath), dataset.location)

    meta = {
        "location": location,
        "slide": {"id": slide.id, "origin_id": slide.origin_id},
        "acquisition": {"id": acquisition.id, "origin_id": acquisition.origin_id},
    }
    return meta, slide_name, acquisition_origin_id


def _import_regionprops(filepath: Path, acquisition_id_mapping: Dict[str, int]):
    p = re.compile("(?P<Name>.*)_(?P<AcquisitionID>[0-9]+).csv")
    slide_name, acquisition_origin_id = p.findall(filepath.name)[0]

    df = pd.read_csv(filepath)
    df.rename(columns={"Object": "ObjectNumber", "centroid-0": "CentroidY", "centroid-1": "CentroidX"}, inplace=True)
    df["ImageNumber"] = acquisition_origin_id
    df["AcquisitionId"] = acquisition_id_mapping.get(f"{slide_name}_{acquisition_origin_id}")
    return slide_name, acquisition_origin_id, df[["ObjectNumber", "ImageNumber", "AcquisitionId", "CentroidX", "CentroidY"]]


def _import_image_csv(src_folder: Path):
    src_uri = src_folder / IMAGE_CSV_FILENAME

    if not src_uri.exists():
        return None

    df = pd.read_csv(src_uri)
    return df


def _import_cell_csv(
    src_folder: Path,
    dst_folder: Path,
    image_number_to_acquisition_id: Dict[int, int],
    channel_order: Dict[str, int],
):
    src_uri = src_folder / CELL_CSV_FILENAME

    if not src_uri.exists():
        return None

    df = pd.read_csv(src_uri)
    df.index = df.index.astype(str, copy=False)

    obs = pd.DataFrame(index=df.index)
    obs["CellId"] = df.index
    obs["AcquisitionId"] = df["ImageNumber"]
    obs["AcquisitionId"].replace(image_number_to_acquisition_id, inplace=True)
    obs["ImageNumber"] = df["ImageNumber"]
    obs["ObjectNumber"] = df["ObjectNumber"]
    obs["CentroidX"] = df["Location_Center_X"]
    obs["CentroidY"] = df["Location_Center_Y"]

    # TODO: skip neighbors columns to keep things simple
    # neighbors_cols = [col for col in df.columns if "Neighbors_" in col]
    # for col in neighbors_cols:
    #     col_name = col.split("_")[1]
    #     if col_name:
    #         obs[col_name] = df[col]
    #         if col_name == "NumberOfNeighbors":
    #             obs[col_name] = obs[col_name].astype("int64")

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
