import logging
import os

from imctools.converters import v1_to_v2
from imctools.io.utils import SESSION_JSON_SUFFIX
from sqlalchemy.orm import Session

from histocat.core.utils import timeit
from histocat.worker.io.imcfolder import import_imcfolder
from histocat.worker.io.utils import locate

logger = logging.getLogger(__name__)


@timeit
def import_imcfolder_v1(db: Session, uri: str, project_id: int):
    input_dir = os.path.dirname(uri)
    output_dir = os.path.join(input_dir, "output")
    v1_to_v2(input_dir, output_dir)

    for session_filename in locate(output_dir, f"*{SESSION_JSON_SUFFIX}"):
        import_imcfolder(db, session_filename, project_id)
