import logging
import os
from typing import List, Optional

from fastapi.encoders import jsonable_encoder
from sqlalchemy.orm import Session

from .dto import SlideCreateDto
from .models import SLIDE_LOCATION_FORMAT, SlideModel

logger = logging.getLogger(__name__)


def get(session: Session, *, id: int) -> Optional[SlideModel]:
    return session.query(SlideModel).filter(SlideModel.id == id).first()


def get_by_name(session: Session, *, experiment_id: int, name: str) -> Optional[SlideModel]:
    return session.query(SlideModel).filter(SlideModel.experiment_id == experiment_id, SlideModel.name == name).first()


def get_multi(session: Session, *, skip: int = 0, limit: int = 100) -> List[SlideModel]:
    return session.query(SlideModel).offset(skip).limit(limit).all()


def create(session: Session, *, params: SlideCreateDto) -> SlideModel:
    data = jsonable_encoder(params)
    entity = SlideModel(**data)
    session.add(entity)
    session.commit()
    session.refresh(entity)

    entity.location = os.path.join(entity.experiment.slides_location, SLIDE_LOCATION_FORMAT.format(id=entity.id))
    if not os.path.exists(entity.location):
        logger.debug(f"Create location for slide {entity.id}: {entity.location}")
        os.makedirs(entity.location)

    session.commit()
    session.refresh(entity)
    return entity


def remove(session: Session, *, id: int):
    item = session.query(SlideModel).filter(SlideModel.id == id).first()
    session.delete(item)
    session.commit()
    return item
