import logging
import zipfile
from pathlib import Path

from imctools.io.utils import MCD_FILENDING, SCHEMA_XML_SUFFIX
from sqlalchemy.orm import Session

from app.core.utils import timeit
from app.io import mcd
from app.io.dataset import CELL_FILENAME, CSV_FILE_EXTENSION, import_dataset
from app.io.imcfolder import import_imcfolder
from app.io.utils import locate

logger = logging.getLogger(__name__)


@timeit
def import_zip(db: Session, uri: str, experiment_id: int, user_id: int):
    path = Path(uri)
    output_dir = path.parent / "output"
    with zipfile.ZipFile(path, "r") as zip:
        zip.extractall(output_dir)

    for mcd_filename in locate(output_dir, f"*{MCD_FILENDING}"):
        mcd.import_mcd(db, mcd_filename, experiment_id, user_id)

    for schema_filename in locate(output_dir, f"*{SCHEMA_XML_SUFFIX}"):
        import_imcfolder(db, schema_filename, experiment_id, user_id)

    for cell_csv_filename in locate(output_dir, f"{CELL_FILENAME}{CSV_FILE_EXTENSION}"):
        import_dataset(db, output_dir, cell_csv_filename, experiment_id, user_id)
