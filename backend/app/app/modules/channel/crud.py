from typing import List, Optional

from fastapi.encoders import jsonable_encoder
from sqlalchemy.orm import Session

from .db import Channel
from .models import ChannelInCreateModel, ChannelInUpdateModel


def get(db_session: Session, *, id: int) -> Optional[Channel]:
    return db_session.query(Channel).filter(Channel.id == id).first()


def get_by_name(db_session: Session, *, name: str) -> Optional[Channel]:
    return db_session.query(Channel).filter(Channel.name == name).first()


def get_multi(db_session: Session, *, skip: int = 0, limit: int = 100) -> List[Optional[Channel]]:
    return db_session.query(Channel).offset(skip).limit(limit).all()


def create(db_session: Session, *, params: ChannelInCreateModel) -> Channel:
    entity = Channel(
        name=params.name,
        metal=params.metal,
        acquisition_id=params.acquisition_id,
        meta=params.meta
    )
    db_session.add(entity)
    db_session.commit()
    db_session.refresh(entity)
    return entity


def update(db_session: Session, *, item: Channel, params: ChannelInUpdateModel) -> Channel:
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
