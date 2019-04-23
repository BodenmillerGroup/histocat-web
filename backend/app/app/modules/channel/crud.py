from typing import List, Optional

from fastapi.encoders import jsonable_encoder
from sqlalchemy.orm import Session

from .db import Channel
from .models import ChannelCreateModel, ChannelUpdateModel


def get(db_session: Session, *, id: int) -> Optional[Channel]:
    return db_session.query(Channel).filter(Channel.id == id).first()


def get_by_name(db_session: Session, *, name: str) -> Optional[Channel]:
    return db_session.query(Channel).filter(Channel.name == name).first()


def get_multi(db_session: Session, *, skip: int = 0, limit: int = 100) -> List[Optional[Channel]]:
    return db_session.query(Channel).offset(skip).limit(limit).all()


def create(db_session: Session, *, params: ChannelCreateModel) -> Channel:
    data = jsonable_encoder(params)
    entity = Channel(**data)
    db_session.add(entity)
    db_session.commit()
    db_session.refresh(entity)
    return entity


def update(db_session: Session, *, item: Channel, params: ChannelUpdateModel) -> Channel:
    data = jsonable_encoder(item)
    update_data = params.dict(skip_defaults=True)
    for field in data:
        if field in update_data:
            setattr(item, field, update_data[field])
    db_session.add(item)
    db_session.commit()
    db_session.refresh(item)
    return item


def remove(db_session: Session, *, id: int):
    item = db_session.query(Channel).filter(Channel.id == id).first()
    db_session.delete(item)
    db_session.commit()
    return item
