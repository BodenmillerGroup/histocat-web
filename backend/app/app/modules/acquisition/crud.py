from typing import List, Optional

from fastapi.encoders import jsonable_encoder
from sqlalchemy.orm import Session

from .db import Acquisition
from .models import AcquisitionInCreateModel, AcquisitionInUpdateModel


def get(db_session: Session, *, id: int) -> Optional[Acquisition]:
    return db_session.query(Acquisition).filter(Acquisition.id == id).first()


def get_by_name(db_session: Session, *, name: str) -> Optional[Acquisition]:
    return db_session.query(Acquisition).filter(Acquisition.name == name).first()


def get_multi(db_session: Session, *, skip: int = 0, limit: int = 100) -> List[Optional[Acquisition]]:
    return db_session.query(Acquisition).offset(skip).limit(limit).all()


def create(db_session: Session, *, params: AcquisitionInCreateModel) -> Acquisition:
    entity = Acquisition(
        name=params.name,
        slide_id=params.slide_id,
        description=params.description,
        meta=params.meta
    )
    db_session.add(entity)
    db_session.commit()
    db_session.refresh(entity)
    return entity


def update(db_session: Session, *, item: Acquisition, params: AcquisitionInUpdateModel) -> Acquisition:
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
