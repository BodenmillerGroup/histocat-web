import logging
import os

from imctools.converters import mcdfolder_to_imcfolder
from imctools.io.utils import SESSION_JSON_SUFFIX
from sqlalchemy.orm import Session

from app.core.utils import timeit
from app.io.imcfolder import import_imcfolder
from app.io.utils import locate

logger = logging.getLogger(__name__)


@timeit
def import_mcd(db: Session, uri: str, experiment_id: int, user_id: int):
    input_dir = os.path.dirname(uri)
    output_dir = os.path.join(input_dir, "output")
    mcdfolder_to_imcfolder(input_dir, output_dir, create_zip=False, skip_csv=True)

    for session_filename in locate(output_dir, f"*{SESSION_JSON_SUFFIX}"):
        import_imcfolder(db, session_filename, experiment_id, user_id)
