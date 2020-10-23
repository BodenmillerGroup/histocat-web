from __future__ import annotations

import logging
import os
import re
from pathlib import Path
from typing import Dict

import anndata as ad
import numpy as np
import pandas as pd
from sqlalchemy.orm import Session

from histocat.core.notifier import Message
from histocat.core.redis_manager import UPDATES_CHANNEL_NAME, redis_manager
from histocat.io.utils import copy_file, locate
from histocat.modules.acquisition import service as acquisition_service
from histocat.modules.dataset import service as dataset_service
from histocat.modules.dataset.dto import DatasetCreateDto, DatasetUpdateDto
from histocat.modules.dataset.models import DatasetModel
from histocat.modules.project import service as project_service
from histocat.modules.slide import service as slide_service

logger = logging.getLogger(__name__)

CSV_FILE_EXTENSION = ".csv"
FEATHER_FILE_EXTENSION = ".feather"
ANNDATA_FILE_EXTENSION = ".h5ad"

EXPERIMENT_CSV_FILENAME = "Experiment"
ACQUISITION_METADATA_FILENAME = "acquisition_metadata"
OBJECT_RELATIONSHIPS_FILENAME = "Object relationships"
IMAGE_FILENAME = "Image"
CELL_FILENAME = "cell"

CHANNELS_FULL_CSV_ENDING = "_ac_full.csv"
PROBABILITIES_MASK_TIFF_ENDING = "_Probabilities_mask.tiff"


def import_dataset(db: Session, root_folder: Path, cell_csv_filename: str, project_id: int):
    """Import dataset from the folder compatible with 'cpout' IMC pipeline folders."""

    project = project_service.get(db, id=project_id)
    if not project:
        logger.warning(f"Cannot import dataset: project [id: {project_id}] does not exist.")
        return

    create_params = DatasetCreateDto(project_id=project_id, status="importing")
    dataset = dataset_service.create(db, params=create_params)
    meta = {}

    src_folder = Path(cell_csv_filename).parent
    dst_folder = Path(dataset.location)
    os.makedirs(dst_folder, exist_ok=True)

    image_df, image_input = _import_image_csv(db, src_folder, dst_folder)
    if image_input:
        meta["image"] = image_input

    probability_masks = {}
    image_map = {}
    image_number_to_acquisition_id = {}
    image_number_to_scaling = {}
    for index, row in image_df.iterrows():
        mask_meta = _import_probabilities_mask(db, src_folder, row, dataset)
        if mask_meta is not None:
            acquisition_id = mask_meta.get("acquisition").get("id")
            probability_masks[acquisition_id] = mask_meta
            image_number = mask_meta.get("image_number")
            image_map[acquisition_id] = image_number
            image_number_to_acquisition_id[image_number] = acquisition_id

            scaling = int(row["Scaling_FullStack"])
            image_number_to_scaling[image_number] = scaling
    meta["probability_masks"] = probability_masks
    # Map acquisition database ID to ImageNumber column
    meta["image_map"] = image_map

    channels_input = {}
    for channels_filename in locate(str(root_folder), f"*{CHANNELS_FULL_CSV_ENDING}"):
        channels_input = _import_channels_csv(db, channels_filename)
        meta["channel_map"] = channels_input
        break

    object_relationships_df, object_relationships_input = _import_object_relationships(db, src_folder, dst_folder)
    if object_relationships_input:
        meta["object_relationships"] = object_relationships_input

    cell_df, cell_input = _import_cell_csv(
        db, src_folder, dst_folder, image_number_to_acquisition_id, image_number_to_scaling, channels_input
    )
    if cell_input:
        meta["cell"] = cell_input

    # Register heatmap columns
    neighbors_cols = [col.split("_")[1] for col in cell_df.columns if "Neighbors_" in col]
    meta["neighbors_columns"] = neighbors_cols

    acquisition_metadata_df, acquisition_metadata_input = _import_acquisition_metadata_csv(db, src_folder, dst_folder)
    if acquisition_metadata_input:
        meta["acquisition_metadata"] = acquisition_metadata_input

    update_params = DatasetUpdateDto(name=f"Dataset {dataset.id}", status="ready", meta=meta)
    dataset = dataset_service.update(db, item=dataset, params=update_params)
    redis_manager.publish(UPDATES_CHANNEL_NAME, Message(project_id, "dataset_imported"))


def _import_image_csv(db: Session, src_folder: Path, dst_folder: Path):
    src_uri = src_folder / f"{IMAGE_FILENAME}{CSV_FILE_EXTENSION}"

    if not src_uri.exists():
        return None, None

    dst_uri = dst_folder / f"{IMAGE_FILENAME}{FEATHER_FILE_EXTENSION}"
    df = pd.read_csv(src_uri)
    df.to_feather(dst_uri)
    artifact = {"location": str(dst_uri)}
    return df, artifact


def _import_channels_csv(db: Session, src_uri: str):
    df = pd.read_csv(src_uri, header=None)
    artifact = {}
    for index, row in df.iterrows():
        artifact[row[0]] = index + 1
    return artifact


def _import_object_relationships(db: Session, src_folder: Path, dst_folder: Path):
    src_uri = src_folder / f"{OBJECT_RELATIONSHIPS_FILENAME}{CSV_FILE_EXTENSION}"

    if not src_uri.exists():
        return None, None

    dst_uri = dst_folder / f"{OBJECT_RELATIONSHIPS_FILENAME}{FEATHER_FILE_EXTENSION}"
    dtype = {
        "Module": "category",
        "Relationship": "category",
        "First Object Name": "category",
        "Second Object Name": "category",
        "First Image Number": "uint16",
        "Second Image Number": "uint16",
        "ModuleNumber": "uint8",
    }
    df = pd.read_csv(src_uri, dtype=dtype)
    df.columns = df.columns.str.replace(" ", "_")
    df.to_feather(dst_uri)
    artifact = {"location": str(dst_uri)}
    return df, artifact


def _import_cell_csv(
    db: Session,
    src_folder: Path,
    dst_folder: Path,
    image_number_to_acquisition_id: Dict[int, int],
    image_number_to_scaling: Dict[int, int],
    channels_input: Dict[str, int],
):
    src_uri = src_folder / f"{CELL_FILENAME}{CSV_FILE_EXTENSION}"

    if not src_uri.exists():
        return None, None

    dst_uri = dst_folder / f"{CELL_FILENAME}{ANNDATA_FILE_EXTENSION}"
    df = pd.read_csv(src_uri)

    obs = pd.DataFrame()
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
    for key, value in channels_input.items():
        # TODO: check intensity multiplier
        x_df[key] = df[f"Intensity_MeanIntensity_FullStackFiltered_c{value}"] * 2 ** 16
        var_names.append(key)
    var = pd.DataFrame(index=var_names)
    var["Channel"] = var.index
    X_counts = x_df.to_numpy()

    adata = ad.AnnData(X_counts, obs=obs, var=var, dtype="float32")
    adata.write_h5ad(dst_uri)

    artifact = {"location": str(dst_uri)}
    return df, artifact


def _import_probabilities_mask(db: Session, src_folder: Path, row: pd.Series, dataset: DatasetModel):
    filename = row["FileName_CellImage"]
    image_number = row["ImageNumber"]
    uri = src_folder / filename

    p = re.compile(
        "(?P<Name>.*)_s(?P<SlideID>[0-9]+)_p(?P<PanoramaID>[0-9]+)_r(?P<AcquisitionROIID>[0-9]+)_a(?P<AcquisitionID>[0-9]+)_ac.*"
    )
    name, slide_origin_id, panorama_origin_id, roi_origin_id, acquisition_origin_id = p.findall(filename)[0]

    slide = slide_service.get_by_name(db, project_id=dataset.project_id, name=name)
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


def _import_acquisition_metadata_csv(db: Session, src_folder: Path, dst_folder: Path):
    src_uri = src_folder / f"{ACQUISITION_METADATA_FILENAME}{CSV_FILE_EXTENSION}"

    if not src_uri.exists():
        return None, None

    dst_uri = dst_folder / f"{ACQUISITION_METADATA_FILENAME}{FEATHER_FILE_EXTENSION}"
    dtype = {
        "MovementType": "category",
        "SegmentDataFormat": "category",
        "SignalType": "category",
        "ValueBytes": "uint8",
        "OrderNumber": "uint8",
        "AcquisitionID": "uint8",
    }
    df = pd.read_csv(src_uri, dtype=dtype)
    df.to_feather(dst_uri)
    artifact = {"location": str(dst_uri)}
    return df, artifact
