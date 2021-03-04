import logging
import zipfile
from pathlib import Path

from sqlalchemy.orm import Session

from histocat.core.model import service as model_service
from histocat.core.utils import timeit

logger = logging.getLogger(__name__)


@timeit
def import_model_zip(db: Session, uri: str, model_id: int):
    model = model_service.get(db, id=model_id)

    path = Path(uri)
    with zipfile.ZipFile(path, "r") as zip:
        zip.extractall(model.location)
