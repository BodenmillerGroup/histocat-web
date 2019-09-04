from __future__ import annotations

import logging
import os
import re
from pathlib import Path

import pandas as pd
from sqlalchemy.orm import Session

from app.core.notifier import Message
from app.core.redis_manager import redis_manager, UPDATES_CHANNEL_NAME
from app.io.utils import copy_file, locate
from app.modules.acquisition import crud as acquisition_crud
from app.modules.dataset import crud as dataset_crud
from app.modules.dataset.db import Dataset
from app.modules.dataset.models import DatasetCreateModel, DatasetUpdateModel
from app.modules.experiment import crud as experiment_crud
from app.modules.panorama import crud as panorama_crud
from app.modules.roi import crud as roi_crud
from app.modules.slide import crud as slide_crud

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


def import_dataset(db: Session, root_folder: Path, cell_csv_filename: str, experiment_id: int, user_id: int):
    """
    Import artifacts from the folder compatible with 'cpout' IMC pipeline folders
    """

    experiment = experiment_crud.get(db, id=experiment_id)
    if not experiment:
        logger.warn(f"Cannot import dataset: experiment [id: {experiment_id}] does not exist.")
        return

    create_params = DatasetCreateModel(
        experiment_id=experiment_id,
        user_id=user_id,
        status="pending",
    )
    dataset = dataset_crud.create(db, params=create_params)
    artifacts = {}

    src_folder = Path(cell_csv_filename).parent
    dst_folder = Path(dataset.location)
    os.makedirs(dst_folder, exist_ok=True)

    image_df, image_artifact = _import_image_csv(db, src_folder, dst_folder)
    if image_artifact:
        artifacts["image"] = image_artifact

    probability_masks = {}
    image_map = {}
    for index, row in image_df.iterrows():
        mask_meta = _import_probabilities_mask(db, src_folder, row, dataset)
        acquisition_id = mask_meta.get("acquisition").get("id")
        probability_masks[acquisition_id] = mask_meta
        image_map[acquisition_id] = mask_meta.get("image_number")
    artifacts["probability_masks"] = probability_masks
    artifacts["image_map"] = image_map

    for channels_filename in locate(str(root_folder), f"*{CHANNELS_FULL_CSV_ENDING}"):
        channels_artifact = _import_channels_csv(db, channels_filename)
        artifacts["channel_map"] = channels_artifact
        break

    object_relationships_df, object_relationships_artifact = _import_object_relationships(db, src_folder, dst_folder)
    if object_relationships_artifact:
        artifacts["object_relationships"] = object_relationships_artifact

    cell_df, cell_artifact = _import_cell_csv(db, src_folder, dst_folder)
    if cell_artifact:
        artifacts["cell"] = cell_artifact

    acquisition_metadata_df, acquisition_metadata_artifact = _import_acquisition_metadata_csv(db, src_folder, dst_folder)
    if acquisition_metadata_artifact:
        artifacts["acquisition_metadata"] = acquisition_metadata_artifact

    update_params = DatasetUpdateModel(
        status="ready",
        artifacts=artifacts,
    )
    dataset = dataset_crud.update(db, item=dataset, params=update_params)
    redis_manager.publish(UPDATES_CHANNEL_NAME, Message(experiment_id, "dataset_imported"))


def _import_image_csv(db: Session, src_folder: Path, dst_folder: Path):
    src_uri = src_folder / f"{IMAGE_FILENAME}{CSV_FILE_EXTENSION}"

    if not src_uri.exists():
        return None, None

    dst_uri = dst_folder / f"{IMAGE_FILENAME}{FEATHER_FILE_EXTENSION}"
    df = pd.read_csv(src_uri)
    df.to_feather(dst_uri)
    artifact = {
        "location": str(dst_uri),
    }
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
    artifact = {
        "location": str(dst_uri),
    }
    return df, artifact


def _import_cell_csv(db: Session, src_folder: Path, dst_folder: Path):
    src_uri = src_folder / f"{CELL_FILENAME}{CSV_FILE_EXTENSION}"

    if not src_uri.exists():
        return None, None

    dst_uri = dst_folder / f"{CELL_FILENAME}{FEATHER_FILE_EXTENSION}"
    df = pd.read_csv(src_uri)
    df.to_feather(dst_uri)
    artifact = {
        "location": str(dst_uri),
    }
    return df, artifact


def _import_probabilities_mask(db: Session, src_folder: Path, row: pd.Series, dataset: Dataset):
    filename = row["FileName_CellImage"]
    image_number = row["ImageNumber"]
    uri = src_folder / filename

    p = re.compile(
        '(?P<Name>.*)_s(?P<SlideID>[0-9]+)_p(?P<PanoramaID>[0-9]+)_r(?P<AcquisitionROIID>[0-9]+)_a(?P<AcquisitionID>[0-9]+)_ac.*')
    name, slide_origin_id, panorama_origin_id, roi_origin_id, acquisition_origin_id = p.findall(filename)[0]

    slide = slide_crud.get_by_name(db, experiment_id=dataset.experiment_id, name=name)
    panorama = panorama_crud.get_by_origin_id(db, slide_id=slide.id, origin_id=panorama_origin_id)
    roi = roi_crud.get_by_origin_id(db, panorama_id=panorama.id, origin_id=roi_origin_id)
    acquisition = acquisition_crud.get_by_origin_id(db, roi_id=roi.id, origin_id=acquisition_origin_id)

    location = copy_file(uri, dataset.location)

    meta = {
        "image_number": image_number,
        "location": location,
        "slide": {
            "id": slide.id,
            "origin_id": slide.origin_id,
        },
        "panorama": {
            "id": panorama.id,
            "origin_id": panorama.origin_id,
        },
        "roi": {
            "id": roi.id,
            "origin_id": roi.origin_id,
        },
        "acquisition": {
            "id": acquisition.id,
            "origin_id": acquisition.origin_id,
        },
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
    artifact = {
        "location": str(dst_uri),
    }
    return df, artifact

