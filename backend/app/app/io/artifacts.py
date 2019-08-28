from __future__ import annotations

import logging
import os
import re
from pathlib import Path

import pandas as pd
from sqlalchemy.orm import Session

from app.io.utils import locate, copy_file
from app.modules.acquisition import crud as acquisition_crud
from app.modules.acquisition_artifact import crud as acquisition_artifact_crud
from app.modules.acquisition_artifact.models import AcquisitionArtifactCreateModel
from app.modules.panorama import crud as panorama_crud
from app.modules.roi import crud as roi_crud
from app.modules.slide import crud as slide_crud
from app.modules.experiment import crud as experiment_crud

logger = logging.getLogger(__name__)

CSV_FILE_EXTENSION = ".csv"
FEATHER_FILE_EXTENSION = ".feather"

EXPERIMENT_CSV_FILENAME = "Experiment"
ACQUISITION_METADATA_FILENAME = "acquisition_metadata"
OBJECT_RELATIONSHIPS_FILENAME = "Object relationships"
IMAGE_FILENAME = "Image"
CELL_FILENAME = "cell"

PROBABILITIES_MASK_TIFF_ENDING = "_Probabilities_mask.tiff"

ARTIFACTS_FOLDER = "artifacts"


def import_artifacts(db: Session, cell_csv_filename: str, experiment_id: int):
    """
    Import artifacts from the folder compatible with 'cpout' IMC pipeline folders
    """

    experiment = experiment_crud.get(db, id=experiment_id)
    if not experiment:
        logger.warn(f"Cannot import artifacts: experiment [id: {experiment_id}] does not exist.")
        return

    src_folder = Path(cell_csv_filename).parent
    dst_folder = Path(experiment.location) / ARTIFACTS_FOLDER
    os.makedirs(dst_folder, exist_ok=True)

    acquisition_metadata_csv = _import_acquisition_metadata(db, src_folder, dst_folder)
    _import_image_csv(db, src_folder, dst_folder)
    _import_object_relationships(db, src_folder, dst_folder)
    _import_cell_csv(db, src_folder, dst_folder)

    for basename in acquisition_metadata_csv['AcSession'].drop_duplicates():
        for mask_filename in locate(src_folder, f"{basename}*{PROBABILITIES_MASK_TIFF_ENDING}"):
            _import_probabilities_mask(db, mask_filename, experiment_id, basename)


def _import_acquisition_metadata(db: Session, src_folder: Path, dst_folder: Path):
    src_uri = src_folder / f"{ACQUISITION_METADATA_FILENAME}{CSV_FILE_EXTENSION}"
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
    return df


def _import_image_csv(db: Session, src_folder: Path, dst_folder: Path):
    src_uri = src_folder / f"{IMAGE_FILENAME}{CSV_FILE_EXTENSION}"
    dst_uri = dst_folder / f"{IMAGE_FILENAME}{FEATHER_FILE_EXTENSION}"
    df = pd.read_csv(src_uri)
    df.to_feather(dst_uri)
    return df


def _import_object_relationships(db: Session, src_folder: Path, dst_folder: Path):
    src_uri = src_folder / f"{OBJECT_RELATIONSHIPS_FILENAME}{CSV_FILE_EXTENSION}"
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
    return df


def _import_cell_csv(db: Session, src_folder: Path, dst_folder: Path):
    src_uri = src_folder / f"{CELL_FILENAME}{CSV_FILE_EXTENSION}"
    dst_uri = dst_folder / f"{CELL_FILENAME}{FEATHER_FILE_EXTENSION}"
    df = pd.read_csv(src_uri)
    df.to_feather(dst_uri)
    return df


def _import_probabilities_mask(db: Session, mask_file_uri: str, experiment_id: int, basename: str):
    path = Path(mask_file_uri)
    filename = path.stem

    slide_original_id = re.findall(f"{basename}_s(\d+)", filename)[0]
    slide_metaname = f"{basename}_s{slide_original_id}"
    slide = slide_crud.get_by_metaname(db, experiment_id=experiment_id, metaname=slide_metaname)

    panorama_original_id = re.findall(f"{slide_metaname}_p(\d+)", filename)[0]
    panorama_metaname = f"{slide_metaname}_p{panorama_original_id}"
    panorama = panorama_crud.get_by_metaname(db, slide_id=slide.id, metaname=panorama_metaname)

    roi_original_id = re.findall(f"{panorama_metaname}_r(\d+)", filename)[0]
    roi_metaname = f"{panorama_metaname}_r{roi_original_id}"
    roi = roi_crud.get_by_metaname(db, panorama_id=panorama.id, metaname=roi_metaname)

    acquisition_original_id = re.findall(f"{roi_metaname}_a(\d+)", filename)[0]
    acquisition_metaname = f"{roi_metaname}_a{acquisition_original_id}"
    acquisition = acquisition_crud.get_by_metaname(db, roi_id=roi.id, metaname=acquisition_metaname)

    dst_folder = os.path.join(slide.location, ARTIFACTS_FOLDER)
    os.makedirs(dst_folder, exist_ok=True)
    location = copy_file(mask_file_uri, dst_folder)

    meta = {
        "slide": {
            "id": slide.id,
            "original_id": slide.original_id,
        },
        "panorama": {
            "id": panorama.id,
            "original_id": panorama.original_id,
        },
        "roi": {
            "id": roi.id,
            "original_id": roi.original_id,
        },
        "acquisition": {
            "id": acquisition.id,
            "original_id": acquisition.original_id,
        },
    }

    params = AcquisitionArtifactCreateModel(
        acquisition_id=acquisition.id,
        type="probabilities_mask",
        location=location,
        meta=meta,
    )
    acquisition_artifact = acquisition_artifact_crud.create(db, params=params)
