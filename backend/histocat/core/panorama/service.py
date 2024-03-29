from typing import List, Optional

from sqlalchemy.orm import Session

from .dto import PanoramaCreateDto
from .models import PanoramaModel


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
    data = params.dict()
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
