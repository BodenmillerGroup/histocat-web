import logging
import os
from typing import List, Optional

from fastapi.encoders import jsonable_encoder
from sqlalchemy.orm import Session

from .db import Slide, SLIDE_LOCATION_FORMAT
from .models import SlideCreateModel

logger = logging.getLogger(__name__)


def get(session: Session, *, id: int) -> Optional[Slide]:
    return session.query(Slide).filter(Slide.id == id).first()


def get_by_name(session: Session, *, name: str) -> Optional[Slide]:
    return session.query(Slide).filter(Slide.name == name).first()


def get_by_experiment_id(session: Session, *, experimentId: int) -> List[Optional[Slide]]:
    return session.query(Slide).filter(Slide.experiment_id == experimentId).all()


def get_multi(session: Session, *, skip: int = 0, limit: int = 100) -> List[Optional[Slide]]:
    return session.query(Slide).offset(skip).limit(limit).all()


def create(session: Session, *, params: SlideCreateModel) -> Slide:
    data = jsonable_encoder(params)
    entity = Slide(**data)
    session.add(entity)
    session.commit()
    session.refresh(entity)

    entity.location = os.path.join(
        entity.experiment.slides_location,
        SLIDE_LOCATION_FORMAT.format(id=entity.id),
    )
    if not os.path.exists(entity.location):
        logger.debug(f'Create location for slide {entity.id}: {entity.location}')
        os.makedirs(entity.location)

    session.commit()
    session.refresh(entity)

    return entity


def remove(session: Session, *, id: int):
    item = session.query(Slide).filter(Slide.id == id).first()
    session.delete(item)
    session.commit()
    return item
