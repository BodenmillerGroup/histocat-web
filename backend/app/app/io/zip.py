import logging
import zipfile
from pathlib import Path

from sqlalchemy.orm import Session

from app.core.utils import timeit
from app.io.imcfolder import import_imcfolder
from app.io.utils import locate, SCHEMA_XML_ENDING

logger = logging.getLogger(__name__)


@timeit
def import_zip(db: Session, uri: str, experiment_id: int):
    path = Path(uri)
    output_dir = path.parent / 'output'
    with zipfile.ZipFile(path, 'r') as zip:
        zip.extractall(output_dir)

    for schema_filename in locate(output_dir, f"*{SCHEMA_XML_ENDING}"):
        import_imcfolder(db, schema_filename, experiment_id)
