import glob
import logging
import os
import zipfile
from pathlib import Path

from imctools.io.utils import MCD_FILENDING, SCHEMA_XML_SUFFIX, SESSION_JSON_SUFFIX
from sqlalchemy.orm import Session

from histocat.core.utils import timeit
from histocat.worker.io import (
    dataset_masks,
    dataset_v1,
    dataset_v2,
    imcfolder,
    imcfolder_v1,
    mcd,
    ometifffolder,
    steinbock,
    utils,
)
from histocat.worker.io.utils import CELL_CSV_FILENAME

logger = logging.getLogger(__name__)

OME_TIFF_SUFFIX = ".ome.tiff"


@timeit
def import_slide(db: Session, uri: str, project_id: int):
    path = Path(uri)
    output_dir = path.parent / "output"
    with zipfile.ZipFile(path, "r") as zip:
        zip.extractall(output_dir)

    for mcd_filename in utils.locate(output_dir, f"*{MCD_FILENDING}"):
        mcd.import_mcd(db, mcd_filename, project_id)

    session_files = glob.glob(os.path.join(output_dir, f"*{SESSION_JSON_SUFFIX}"))
    schema_files = glob.glob(os.path.join(output_dir, f"*{SCHEMA_XML_SUFFIX}"))
    ome_tiff_files = glob.glob(os.path.join(output_dir, "**", f"*{OME_TIFF_SUFFIX}"), recursive=True)

    if len(session_files) > 0 and len(schema_files) > 0:
        for schema_filename in schema_files:
            imcfolder.import_imcfolder(db, schema_filename, project_id)
    elif len(session_files) == 0 and len(schema_files) > 0:
        for schema_filename in schema_files:
            imcfolder_v1.import_imcfolder_v1(db, schema_filename, project_id)
    elif len(ome_tiff_files) > 0:
        ometifffolder.import_ometifffolder(db, str(path.stem), ome_tiff_files, project_id)


@timeit
def import_dataset(
    db: Session,
    type: str,
    masks_folder: str,
    regionprops_folder: str,
    intensities_folder: str,
    uri: str,
    project_id: int,
):
    path = Path(uri)
    output_dir = path.parent / "output"
    with zipfile.ZipFile(path, "r") as zip:
        zip.extractall(output_dir)

    if type == "steinbock":
        steinbock.import_dataset(db, output_dir, project_id, masks_folder, regionprops_folder, intensities_folder)
    elif type == "ImcSegmentationPipelineV1":
        for cell_csv_filename in utils.locate(output_dir, CELL_CSV_FILENAME):
            dataset_v1.import_dataset(db, output_dir, cell_csv_filename, project_id)
    elif type == "ImcSegmentationPipelineV2":
        for cell_csv_filename in utils.locate(output_dir, CELL_CSV_FILENAME):
            dataset_v2.import_dataset(db, output_dir, cell_csv_filename, project_id)
    elif type == "masks":
        masks_csv_files = glob.glob(os.path.join(output_dir, "**", dataset_masks.MASKS_CSV_FILE), recursive=True)
        if len(masks_csv_files) == 1:
            dataset_masks.import_dataset(db, Path(masks_csv_files[0]), project_id)
