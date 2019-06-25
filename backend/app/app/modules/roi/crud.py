import logging
import os
from typing import List, Optional

from fastapi.encoders import jsonable_encoder
from sqlalchemy.orm import Session

from .db import ROI, ROI_LOCATION_FORMAT
from .models import ROICreateModel

logger = logging.getLogger(__name__)


def get(session: Session, *, id: int) -> Optional[ROI]:
    return session.query(ROI).filter(ROI.id == id).first()


def get_by_panorama_id(session: Session, *, panorama_id: int) -> List[ROI]:
    return session.query(ROI).filter(ROI.panorama_id == ROI).all()


def get_multi(session: Session, *, skip: int = 0, limit: int = 100) -> List[ROI]:
    return session.query(ROI).offset(skip).limit(limit).all()


def create(session: Session, *, params: ROICreateModel) -> ROI:
    data = jsonable_encoder(params)
    entity = ROI(**data)
    session.add(entity)
    session.commit()
    session.refresh(entity)

    entity.location = os.path.join(
        entity.panorama.rois_location,
        ROI_LOCATION_FORMAT.format(id=entity.id),
    )
    if not os.path.exists(entity.location):
        logger.debug(f'Create location for ROI {entity.id}: {entity.location}')
        os.makedirs(entity.location)

    session.commit()
    session.refresh(entity)
    return entity


def remove(session: Session, *, id: int):
    item = session.query(ROI).filter(ROI.id == id).first()
    session.delete(item)
    session.commit()
    return item
