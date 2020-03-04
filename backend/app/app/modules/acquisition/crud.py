import logging
from typing import List, Optional

from fastapi.encoders import jsonable_encoder
from sqlalchemy.orm import Session

from .db import Acquisition
from .models import AcquisitionCreateModel

logger = logging.getLogger(__name__)


def get(session: Session, *, id: int) -> Optional[Acquisition]:
    return session.query(Acquisition).filter(Acquisition.id == id).first()


def get_by_origin_id(session: Session, *, roi_id: int, origin_id: int) -> Optional[Acquisition]:
    return session.query(Acquisition).filter(Acquisition.roi_id == roi_id, Acquisition.origin_id == origin_id).first()


def get_multi(session: Session, *, skip: int = 0, limit: int = 1000) -> List[Optional[Acquisition]]:
    return session.query(Acquisition).offset(skip).limit(limit).all()


def create(session: Session, *, params: AcquisitionCreateModel) -> Acquisition:
    data = jsonable_encoder(params)
    entity = Acquisition(**data)
    session.add(entity)
    session.commit()
    session.refresh(entity)
    return entity


def remove(session: Session, *, id: int):
    item = session.query(Acquisition).filter(Acquisition.id == id).first()
    session.delete(item)
    session.commit()
    return item
