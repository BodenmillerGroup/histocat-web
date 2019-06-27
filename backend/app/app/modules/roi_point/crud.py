import logging
from typing import List, Optional

from fastapi.encoders import jsonable_encoder
from sqlalchemy.orm import Session

from .db import ROIPoint
from .models import ROIPointCreateModel

logger = logging.getLogger(__name__)


def get(session: Session, *, id: int) -> Optional[ROIPoint]:
    return session.query(ROIPoint).filter(ROIPoint.id == id).first()


def get_by_roi_id(session: Session, *, roi_id: int) -> List[ROIPoint]:
    return session.query(ROIPoint).filter(ROIPoint.roi_id == roi_id).all()


def get_multi(session: Session, *, skip: int = 0, limit: int = 100) -> List[ROIPoint]:
    return session.query(ROIPoint).offset(skip).limit(limit).all()


def create(session: Session, *, params: ROIPointCreateModel) -> ROIPoint:
    data = jsonable_encoder(params)
    entity = ROIPoint(**data)
    session.add(entity)
    session.commit()
    session.refresh(entity)
    return entity


def remove(session: Session, *, id: int):
    item = session.query(ROIPoint).filter(ROIPoint.id == id).first()
    session.delete(item)
    session.commit()
    return item
