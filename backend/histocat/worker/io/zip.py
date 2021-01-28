import glob
import logging
import os
import zipfile
from pathlib import Path

from imctools.io.utils import MCD_FILENDING, SCHEMA_XML_SUFFIX, SESSION_JSON_SUFFIX
from sqlalchemy.orm import Session

from histocat.core.constants import CSV_FILE_EXTENSION
from histocat.core.utils import timeit
from histocat.worker.io import mcd
from histocat.worker.io.dataset_v1 import CELL_FILENAME, import_dataset
from histocat.worker.io.imcfolder import import_imcfolder
from histocat.worker.io.imcfolder_v1 import import_imcfolder_v1
from histocat.worker.io.utils import locate

logger = logging.getLogger(__name__)


@timeit
def import_slide_zip(db: Session, uri: str, project_id: int):
    path = Path(uri)
    output_dir = path.parent / "output"
    with zipfile.ZipFile(path, "r") as zip:
        zip.extractall(output_dir)

    for mcd_filename in locate(output_dir, f"*{MCD_FILENDING}"):
        mcd.import_mcd(db, mcd_filename, project_id)

    session_files = glob.glob(os.path.join(output_dir, f"*{SESSION_JSON_SUFFIX}"))
    schema_files = glob.glob(os.path.join(output_dir, f"*{SCHEMA_XML_SUFFIX}"))

    if len(session_files) > 0 and len(schema_files) > 0:
        for schema_filename in schema_files:
            import_imcfolder(db, schema_filename, project_id)
    elif len(session_files) == 0 and len(schema_files) > 0:
        for schema_filename in schema_files:
            import_imcfolder_v1(db, schema_filename, project_id)


@timeit
def import_dataset_zip(db: Session, uri: str, project_id: int):
    path = Path(uri)
    output_dir = path.parent / "output"
    with zipfile.ZipFile(path, "r") as zip:
        zip.extractall(output_dir)

    for cell_csv_filename in locate(output_dir, f"{CELL_FILENAME}{CSV_FILE_EXTENSION}"):
        import_dataset(db, output_dir, cell_csv_filename, project_id)
