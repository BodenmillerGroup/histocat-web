import logging
from typing import Optional

from fastapi.encoders import jsonable_encoder
from sqlalchemy.orm import Session

from .models import Acquisition
from .dto import AcquisitionCreateDto

logger = logging.getLogger(__name__)


def get_by_id(session: Session, id: int) -> Optional[Acquisition]:
    return session.query(Acquisition).filter(Acquisition.id == id).first()


def get_by_origin_id(session: Session, *, roi_id: int, origin_id: int) -> Optional[Acquisition]:
    return session.query(Acquisition).filter(Acquisition.roi_id == roi_id, Acquisition.origin_id == origin_id).first()


def create(session: Session, params: AcquisitionCreateDto) -> Acquisition:
    data = jsonable_encoder(params)
    item = Acquisition(**data)
    session.add(item)
    session.commit()
    session.refresh(item)
    return item


def delete_by_id(session: Session, id: int):
    item = get_by_id(session, id)
    session.delete(item)
    session.commit()
    return item
