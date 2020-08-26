from __future__ import annotations

import logging
import os
import re
from pathlib import Path
from typing import Dict

import pandas as pd
from sqlalchemy.orm import Session

from histocat.core.notifier import Message
from histocat.core.redis_manager import UPDATES_CHANNEL_NAME, redis_manager
from histocat.io.utils import copy_file, locate
from histocat.modules.acquisition import service as acquisition_service
from histocat.modules.dataset import service as dataset_service
from histocat.modules.dataset.dto import DatasetCreateDto, DatasetUpdateDto
from histocat.modules.dataset.models import DatasetModel
from histocat.modules.experiment import service as experiment_service
from histocat.modules.slide import service as slide_service

logger = logging.getLogger(__name__)

CSV_FILE_EXTENSION = ".csv"
FEATHER_FILE_EXTENSION = ".feather"

EXPERIMENT_CSV_FILENAME = "Experiment"
ACQUISITION_METADATA_FILENAME = "acquisition_metadata"
OBJECT_RELATIONSHIPS_FILENAME = "Object relationships"
IMAGE_FILENAME = "Image"
CELL_FILENAME = "cell"

CHANNELS_FULL_CSV_ENDING = "_ac_full.csv"
PROBABILITIES_MASK_TIFF_ENDING = "_Probabilities_mask.tiff"


def import_dataset(db: Session, root_folder: Path, cell_csv_filename: str, experiment_id: int):
    """Import dataset from the folder compatible with 'cpout' IMC pipeline folders."""

    experiment = experiment_service.get(db, id=experiment_id)
    if not experiment:
        logger.warning(f"Cannot import dataset: experiment [id: {experiment_id}] does not exist.")
        return

    create_params = DatasetCreateDto(experiment_id=experiment_id, status="pending")
    dataset = dataset_service.create(db, params=create_params)
    input = {}

    src_folder = Path(cell_csv_filename).parent
    dst_folder = Path(dataset.location)
    os.makedirs(dst_folder, exist_ok=True)

    image_df, image_input = _import_image_csv(db, src_folder, dst_folder)
    if image_input:
        input["image"] = image_input

    probability_masks = {}
    image_map = {}
    image_number_to_acquisition_id = {}
    for index, row in image_df.iterrows():
        mask_meta = _import_probabilities_mask(db, src_folder, row, dataset)
        if mask_meta is not None:
            acquisition_id = mask_meta.get("acquisition").get("id")
            probability_masks[acquisition_id] = mask_meta
            image_number = mask_meta.get("image_number")
            image_map[acquisition_id] = image_number
            image_number_to_acquisition_id[image_number] = acquisition_id
    input["probability_masks"] = probability_masks
    # Map acquisition database ID to ImageNumber column
    input["image_map"] = image_map

    channels_input = {}
    for channels_filename in locate(str(root_folder), f"*{CHANNELS_FULL_CSV_ENDING}"):
        channels_input = _import_channels_csv(db, channels_filename)
        input["channel_map"] = channels_input
        break

    object_relationships_df, object_relationships_input = _import_object_relationships(db, src_folder, dst_folder)
    if object_relationships_input:
        input["object_relationships"] = object_relationships_input

    cell_df, cell_input = _import_cell_csv(db, src_folder, dst_folder, image_number_to_acquisition_id, channels_input)
    if cell_input:
        input["cell"] = cell_input

    # Register heatmap columns
    neighbors_cols = [col.split("_")[1] for col in cell_df.columns if "Neighbors_" in col]
    input["neighbors_columns"] = neighbors_cols

    acquisition_metadata_df, acquisition_metadata_input = _import_acquisition_metadata_csv(db, src_folder, dst_folder)
    if acquisition_metadata_input:
        input["acquisition_metadata"] = acquisition_metadata_input

    update_params = DatasetUpdateDto(status="ready", input=input)
    dataset = dataset_service.update(db, item=dataset, params=update_params)
    redis_manager.publish(UPDATES_CHANNEL_NAME, Message(experiment_id, "dataset_imported"))


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
    channels_input: Dict[str, int],
):
    src_uri = src_folder / f"{CELL_FILENAME}{CSV_FILE_EXTENSION}"

    if not src_uri.exists():
        return None, None

    dst_uri = dst_folder / f"{CELL_FILENAME}{FEATHER_FILE_EXTENSION}"
    df = pd.read_csv(src_uri)

    df_preprocessed = pd.DataFrame()
    df_preprocessed["CellId"] = df.index
    df_preprocessed["AcquisitionId"] = df["ImageNumber"]
    df_preprocessed["AcquisitionId"].replace(image_number_to_acquisition_id, inplace=True)
    df_preprocessed["ImageNumber"] = df["ImageNumber"]
    df_preprocessed["ObjectNumber"] = df["ObjectNumber"]
    df_preprocessed["CentroidX"] = df["Location_Center_X"]
    df_preprocessed["CentroidY"] = df["Location_Center_Y"]
    for key, value in channels_input.items():
        # TODO: check intensity multiplier
        df_preprocessed[key] = df[f"Intensity_MeanIntensity_FullStack_c{value}"] * 2 ** 16

    neighbors_cols = [col for col in df.columns if "Neighbors_" in col]
    for col in neighbors_cols:
        col_name = col.split("_")[1]
        if col_name:
            df_preprocessed[col_name] = df[col]
            if col_name == "NumberOfNeighbors":
                df_preprocessed[col_name] = df_preprocessed[col_name].astype("int64")

    df_preprocessed.to_feather(dst_uri)
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

    slide = slide_service.get_by_name(db, experiment_id=dataset.experiment_id, name=name)
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
