import logging
from typing import List, Optional

from fastapi.encoders import jsonable_encoder
from sqlalchemy.orm import Session

from .models import PanoramaModel
from .dto import PanoramaCreateDto

logger = logging.getLogger(__name__)


def get(session: Session, *, id: int) -> Optional[PanoramaModel]:
    return session.query(PanoramaModel).filter(PanoramaModel.id == id).first()


def get_by_origin_id(session: Session, *, slide_id: int, origin_id: int) -> Optional[PanoramaModel]:
    return (
        session.query(PanoramaModel)
        .filter(PanoramaModel.slide_id == slide_id, PanoramaModel.origin_id == origin_id)
        .first()
    )


def get_by_slide_id(session: Session, *, slide_id: int) -> List[PanoramaModel]:
    return session.query(PanoramaModel).filter(PanoramaModel.slide_id == slide_id).all()


def get_multi(session: Session, *, skip: int = 0, limit: int = 100) -> List[PanoramaModel]:
    return session.query(PanoramaModel).offset(skip).limit(limit).all()


def create(session: Session, *, params: PanoramaCreateDto) -> PanoramaModel:
    data = jsonable_encoder(params)
    entity = PanoramaModel(**data)
    session.add(entity)
    session.commit()
    session.refresh(entity)
    return entity


def remove(session: Session, *, id: int):
    item = session.query(PanoramaModel).filter(PanoramaModel.id == id).first()
    session.delete(item)
    session.commit()
    return item
