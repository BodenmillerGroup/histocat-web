from typing import List, Optional

from fastapi.encoders import jsonable_encoder
from sqlalchemy.orm import Session

from .db import Acquisition
from .models import AcquisitionCreateModel, AcquisitionUpdateModel


def get(session: Session, *, id: int) -> Optional[Acquisition]:
    return session.query(Acquisition).filter(Acquisition.id == id).first()


def get_by_name(session: Session, *, name: str) -> Optional[Acquisition]:
    return session.query(Acquisition).filter(Acquisition.name == name).first()


def get_multi(session: Session, *, skip: int = 0, limit: int = 100) -> List[Optional[Acquisition]]:
    return session.query(Acquisition).offset(skip).limit(limit).all()


def create(session: Session, *, params: AcquisitionCreateModel) -> Acquisition:
    data = jsonable_encoder(params)
    entity = Acquisition(**data)
    session.add(entity)
    session.commit()
    session.refresh(entity)
    return entity


def update(session: Session, *, item: Acquisition, params: AcquisitionUpdateModel) -> Acquisition:
    data = jsonable_encoder(item)
    update_data = params.dict(skip_defaults=True)
    for field in data:
        if field in update_data:
            setattr(item, field, update_data[field])
    session.add(item)
    session.commit()
    session.refresh(item)
    return item


def remove(session: Session, *, id: int):
    item = session.query(Acquisition).filter(Acquisition.id == id).first()
    session.delete(item)
    session.commit()
    return item
