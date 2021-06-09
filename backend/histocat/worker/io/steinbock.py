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
from histocat.worker.io.utils import copy_file

logger = logging.getLogger(__name__)


PANEL_CSV_FILE = "panel.csv"


def _report_error(project_id: int, message: str):
    """Log error message and send it to the client via websocket"""
    logger.warning(message)
    redis_manager.publish(UPDATES_CHANNEL_NAME, Message(project_id, "error", message))


def import_dataset(db: Session, input_folder: Path, project_id: int, masks_folder: str, regionprops_folder: str, intensities_folder: str):
    """Import dataset from the folder compatible with 'steinbock' format."""

    # Validate data

    project = project_service.get(db, id=project_id)
    if not project:
        _report_error(project_id, f"Dataset Import Error: project [id: {project_id}] does not exist")
        return

    # Find data source folder where panel file resides
    src_folder = None
    for path in input_folder.rglob(PANEL_CSV_FILE):
        src_folder = path.parent
        break

    if src_folder is None:
        _report_error(project_id, f"Dataset Import Error: panel file is missing")
        return

    mask_files = sorted(Path(src_folder / masks_folder).rglob("*.tiff"))
    if len(mask_files) == 0:
        _report_error(project_id, f"Dataset Import Error: mask files are missing in folder '{masks_folder}'")
        return

    regionprops_files = sorted(Path(src_folder / regionprops_folder).rglob("*.csv"))
    if len(regionprops_files) == 0:
        _report_error(project_id, f"Dataset Import Error: regionprops files are missing in folder '{regionprops_folder}'")
        return

    intensities_files = sorted(Path(src_folder / intensities_folder).rglob("*.csv"))
    if len(intensities_files) == 0:
        _report_error(project_id,
                      f"Dataset Import Error: intensities files are missing in folder '{intensities_folder}'")
        return

    # Postpone dataset db entry creation until input data validated
    create_params = DatasetCreateDto(project_id=project_id, origin="steinbock", status="pending")
    dataset = dataset_service.create(db, params=create_params)

    dst_folder = Path(dataset.location)
    os.makedirs(dst_folder, exist_ok=True)

    # Import panel data
    panel_df = _import_panel(os.path.join(src_folder, PANEL_CSV_FILE))

    # Metadata dictionary
    meta = {}
    masks = {}
    acquisition_id_mapping = {}

    for mask_file in mask_files:
        result = _import_mask(db, mask_file, dataset)
        if result is not None:
            mask_meta, slide_name, acquisition_origin_id = result
            acquisition_id = mask_meta.get("acquisition").get("id")
            image_number = mask_meta.get("acquisition").get("origin_id")
            masks[acquisition_id] = mask_meta
            acquisition_id_mapping[f"{slide_name}_{image_number}"] = acquisition_id
    meta["masks"] = masks

    regionprops_df = pd.DataFrame()
    for regionprops_file in regionprops_files:
        slide_name, acquisition_origin_id, df = _import_regionprops(regionprops_file, acquisition_id_mapping)
        regionprops_df = regionprops_df.append(df)
    regionprops_df.reset_index(inplace=True, drop=True)
    regionprops_df["CellId"] = regionprops_df.index
    print(regionprops_df)

    intensities_df = pd.DataFrame()
    for intensities_file in intensities_files:
        slide_name, acquisition_origin_id, df = _import_intensities(intensities_file, acquisition_id_mapping)
        intensities_df = intensities_df.append(df)
    intensities_df.reset_index(inplace=True, drop=True)
    intensities_df["CellId"] = regionprops_df.index
    print(intensities_df)

    var_names = []
    x_df = pd.DataFrame()
    for index, row in panel_df.iterrows():
        channel = row["channel"]
        name = row["name"]
        # TODO: check intensity
        x_df[channel] = intensities_df[name]
        var_names.append(channel)

    var = pd.DataFrame(index=var_names)
    var["Channel"] = var.index
    X_counts = x_df.to_numpy()

    adata = ad.AnnData(X_counts, obs=regionprops_df, var=var, dtype="float32")
    dst_uri = dst_folder / CELL_FILENAME
    adata.write_h5ad(dst_uri)

    acquisition_ids = sorted(list(masks.keys()))

    update_params = DatasetUpdateDto(
        name=f"Dataset {dataset.id}", status="ready", acquisition_ids=acquisition_ids, channels=var_names, meta=meta
    )
    dataset = dataset_service.update(db, item=dataset, params=update_params)
    redis_manager.publish(UPDATES_CHANNEL_NAME, Message(project_id, "dataset_imported"))


def _import_panel(path: str):
    panel_df = pd.read_csv(path)
    return panel_df


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


def _import_intensities(filepath: Path, acquisition_id_mapping: Dict[str, int]):
    p = re.compile("(?P<Name>.*)_(?P<AcquisitionID>[0-9]+).csv")
    slide_name, acquisition_origin_id = p.findall(filepath.name)[0]

    df = pd.read_csv(filepath)
    df.rename(columns={"Object": "ObjectNumber"}, inplace=True)
    df["ImageNumber"] = acquisition_origin_id
    df["AcquisitionId"] = acquisition_id_mapping.get(f"{slide_name}_{acquisition_origin_id}")
    return slide_name, acquisition_origin_id, df
