import logging
import os
from typing import List, Optional

from fastapi.encoders import jsonable_encoder
from sqlalchemy.orm import Session

from .db import Panorama, PANORAMA_LOCATION_FORMAT
from .models import PanoramaCreateModel

logger = logging.getLogger(__name__)


def get(session: Session, *, id: int) -> Optional[Panorama]:
    return session.query(Panorama).filter(Panorama.id == id).first()


def get_by_slide_id(session: Session, *, slide_id: int) -> List[Panorama]:
    return session.query(Panorama).filter(Panorama.slide_id == slide_id).all()


def get_multi(session: Session, *, skip: int = 0, limit: int = 100) -> List[Panorama]:
    return session.query(Panorama).offset(skip).limit(limit).all()


def create(session: Session, *, params: PanoramaCreateModel) -> Panorama:
    data = jsonable_encoder(params)
    entity = Panorama(**data)
    session.add(entity)
    session.commit()
    session.refresh(entity)

    entity.location = os.path.join(
        entity.slide.panoramas_location,
        PANORAMA_LOCATION_FORMAT.format(id=entity.id),
    )
    if not os.path.exists(entity.location):
        logger.debug(f'Create location for panorama {entity.id}: {entity.location}')
        os.makedirs(entity.location)

    session.commit()
    session.refresh(entity)
    return entity


def remove(session: Session, *, id: int):
    item = session.query(Panorama).filter(Panorama.id == id).first()
    session.delete(item)
    session.commit()
    return item
