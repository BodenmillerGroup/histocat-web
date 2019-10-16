import logging
import os

from imctools.scripts.convertfolder2imcfolder import convert_folder2imcfolder
from sqlalchemy.orm import Session

from app.core.utils import timeit
from app.io.imcfolder import import_imcfolder
from app.io.utils import locate, SCHEMA_XML_ENDING

logger = logging.getLogger(__name__)


@timeit
def import_mcd(db: Session, uri: str, experiment_id: int, user_id: int):
    input_dir = os.path.dirname(uri)
    output_dir = os.path.join(input_dir, "output")
    convert_folder2imcfolder(input_dir, output_dir, dozip=False)

    for schema_filename in locate(output_dir, f"*{SCHEMA_XML_ENDING}"):
        import_imcfolder(db, schema_filename, experiment_id, user_id)
