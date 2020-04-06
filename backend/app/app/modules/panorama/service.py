import logging
from typing import List, Optional

from fastapi.encoders import jsonable_encoder
from sqlalchemy.orm import Session

from .models import Panorama
from .dto import PanoramaCreateDto

logger = logging.getLogger(__name__)


def get(session: Session, *, id: int) -> Optional[Panorama]:
    return session.query(Panorama).filter(Panorama.id == id).first()


def get_by_origin_id(session: Session, *, slide_id: int, origin_id: int) -> Optional[Panorama]:
    return session.query(Panorama).filter(Panorama.slide_id == slide_id, Panorama.origin_id == origin_id).first()


def get_by_slide_id(session: Session, *, slide_id: int) -> List[Panorama]:
    return session.query(Panorama).filter(Panorama.slide_id == slide_id).all()


def get_multi(session: Session, *, skip: int = 0, limit: int = 100) -> List[Panorama]:
    return session.query(Panorama).offset(skip).limit(limit).all()


def create(session: Session, *, params: PanoramaCreateDto) -> Panorama:
    data = jsonable_encoder(params)
    entity = Panorama(**data)
    session.add(entity)
    session.commit()
    session.refresh(entity)
    return entity


def remove(session: Session, *, id: int):
    item = session.query(Panorama).filter(Panorama.id == id).first()
    session.delete(item)
    session.commit()
    return item
