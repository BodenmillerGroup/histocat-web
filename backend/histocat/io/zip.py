import glob
import logging
import os
import zipfile
from pathlib import Path

from imctools.io.utils import MCD_FILENDING, SCHEMA_XML_SUFFIX, SESSION_JSON_SUFFIX
from sqlalchemy.orm import Session

from histocat.core.utils import timeit
from histocat.io import mcd
from histocat.io.dataset import CELL_FILENAME, CSV_FILE_EXTENSION, import_dataset
from histocat.io.imcfolder import import_imcfolder
from histocat.io.imcfolder_v1 import import_imcfolder_v1
from histocat.io.utils import locate

logger = logging.getLogger(__name__)


@timeit
def import_zip(db: Session, uri: str, experiment_id: int):
    path = Path(uri)
    output_dir = path.parent / "output"
    with zipfile.ZipFile(path, "r") as zip:
        zip.extractall(output_dir)

    for mcd_filename in locate(output_dir, f"*{MCD_FILENDING}"):
        mcd.import_mcd(db, mcd_filename, experiment_id)

    session_files = glob.glob(os.path.join(output_dir, f"*{SESSION_JSON_SUFFIX}"))
    schema_files = glob.glob(os.path.join(output_dir, f"*{SCHEMA_XML_SUFFIX}"))

    if len(session_files) > 0 and len(schema_files) > 0:
        for schema_filename in schema_files:
            import_imcfolder(db, schema_filename, experiment_id)
    elif len(session_files) == 0 and len(schema_files) > 0:
        for schema_filename in schema_files:
            import_imcfolder_v1(db, schema_filename, experiment_id)

    for cell_csv_filename in locate(output_dir, f"{CELL_FILENAME}{CSV_FILE_EXTENSION}"):
        import_dataset(db, output_dir, cell_csv_filename, experiment_id)
