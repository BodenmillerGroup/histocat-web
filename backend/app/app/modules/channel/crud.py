import logging
import os
from typing import List, Optional

from fastapi.encoders import jsonable_encoder
from sqlalchemy.orm import Session

from .db import Channel, CHANNEL_LOCATION_FORMAT
from .models import ChannelCreateModel

logger = logging.getLogger(__name__)


def get(session: Session, *, id: int) -> Optional[Channel]:
    return session.query(Channel).filter(Channel.id == id).first()


def get_by_name(session: Session, *, name: str) -> Optional[Channel]:
    return session.query(Channel).filter(Channel.name == name).first()


def get_multi(session: Session, *, skip: int = 0, limit: int = 100) -> List[Optional[Channel]]:
    return session.query(Channel).offset(skip).limit(limit).all()


def create(session: Session, *, params: ChannelCreateModel) -> Channel:
    data = jsonable_encoder(params)
    entity = Channel(**data)
    session.add(entity)
    session.commit()
    session.refresh(entity)

    entity.location = os.path.join(
        entity.acquisition.channels_location,
        CHANNEL_LOCATION_FORMAT.format(id=entity.id),
    )
    if not os.path.exists(entity.location):
        logger.debug(f'Create location for channel {entity.id}: {entity.location}')
        os.makedirs(entity.location)

    session.commit()
    session.refresh(entity)

    return entity


def remove(session: Session, *, id: int):
    item = session.query(Channel).filter(Channel.id == id).first()
    session.delete(item)
    session.commit()
    return item
