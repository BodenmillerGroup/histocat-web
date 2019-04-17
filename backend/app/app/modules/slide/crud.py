from typing import List, Optional

from fastapi.encoders import jsonable_encoder
from sqlalchemy.orm import Session

from .db import Slide
from .models import SlideInCreateModel, SlideInUpdateModel


def get(db_session: Session, *, id: int) -> Optional[Slide]:
    return db_session.query(Slide).filter(Slide.id == id).first()


def get_by_name(db_session: Session, *, name: str) -> Optional[Slide]:
    return db_session.query(Slide).filter(Slide.name == name).first()


def get_multi(db_session: Session, *, skip: int = 0, limit: int = 100) -> List[Optional[Slide]]:
    return db_session.query(Slide).offset(skip).limit(limit).all()


def create(db_session: Session, *, params: SlideInCreateModel) -> Slide:
    entity = Slide(
        name=params.name,
        experiment_id=params.experiment_id,
        description=params.description,
        meta=params.meta
    )
    db_session.add(entity)
    db_session.commit()
    db_session.refresh(entity)
    return entity


def update(db_session: Session, *, item: Slide, params: SlideInUpdateModel) -> Slide:
    data = jsonable_encoder(item)
    for field in data:
        if field in params.fields:
            value_in = getattr(params, field)
            if value_in is not None:
                setattr(item, field, value_in)
    db_session.add(item)
    db_session.commit()
    db_session.refresh(item)
    return item
