from typing import List, Optional

from fastapi.encoders import jsonable_encoder
from sqlalchemy.orm import Session

from .db import Channel
from .models import ChannelCreateModel, ChannelUpdateModel


def get(session: Session, *, id: int) -> Optional[Channel]:
    return session.query(Channel).filter(Channel.id == id).first()


def get_by_name(session: Session, *, name: str) -> Optional[Channel]:
    return session.query(Channel).filter(Channel.name == name).first()


def get_multi(
    session: Session, *, skip: int = 0, limit: int = 100
) -> List[Optional[Channel]]:
    return session.query(Channel).offset(skip).limit(limit).all()


def create(session: Session, *, params: ChannelCreateModel) -> Channel:
    data = jsonable_encoder(params)
    entity = Channel(**data)
    session.add(entity)
    session.commit()
    session.refresh(entity)
    return entity


def update(session: Session, *, item: Channel, params: ChannelUpdateModel) -> Channel:
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
    item = session.query(Channel).filter(Channel.id == id).first()
    session.delete(item)
    session.commit()
    return item
