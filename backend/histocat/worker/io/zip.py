import glob
import logging
import os
import zipfile
from pathlib import Path

from sqlalchemy.orm import Session
from imctools.io.utils import MCD_FILENDING, SESSION_JSON_SUFFIX, SCHEMA_XML_SUFFIX

from histocat.core.constants import CSV_FILE_EXTENSION
from histocat.core.utils import timeit
from histocat.worker.io import imcfolder_v1, dataset_v1, dataset_v2, mcd, imcfolder, utils

logger = logging.getLogger(__name__)


@timeit
def import_slide_zip(db: Session, uri: str, project_id: int):
    path = Path(uri)
    output_dir = path.parent / "output"
    with zipfile.ZipFile(path, "r") as zip:
        zip.extractall(output_dir)

    for mcd_filename in utils.locate(output_dir, f"*{MCD_FILENDING}"):
        mcd.import_mcd(db, mcd_filename, project_id)

    session_files = glob.glob(os.path.join(output_dir, f"*{SESSION_JSON_SUFFIX}"))
    schema_files = glob.glob(os.path.join(output_dir, f"*{SCHEMA_XML_SUFFIX}"))

    if len(session_files) > 0 and len(schema_files) > 0:
        for schema_filename in schema_files:
            imcfolder.import_imcfolder(db, schema_filename, project_id)
    elif len(session_files) == 0 and len(schema_files) > 0:
        for schema_filename in schema_files:
            imcfolder_v1.import_imcfolder_v1(db, schema_filename, project_id)


@timeit
def import_dataset_zip(db: Session, uri: str, project_id: int):
    path = Path(uri)
    output_dir = path.parent / "output"
    with zipfile.ZipFile(path, "r") as zip:
        zip.extractall(output_dir)

    for cell_csv_filename in utils.locate(output_dir, f"{dataset_v1.CELL_FILENAME}{CSV_FILE_EXTENSION}"):
        src_folder = Path(cell_csv_filename).parent
        is_v2 = os.path.exists(os.path.join(src_folder, "var_cell.csv"))
        if is_v2:
            dataset_v2.import_dataset(db, output_dir, cell_csv_filename, project_id)
        else:
            dataset_v1.import_dataset(db, output_dir, cell_csv_filename, project_id)
