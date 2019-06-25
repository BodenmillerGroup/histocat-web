import logging
import os
from typing import List, Optional

from fastapi.encoders import jsonable_encoder
from sqlalchemy.orm import Session

from .db import Acquisition, ACQUISITION_LOCATION_FORMAT
from .models import AcquisitionCreateModel

logger = logging.getLogger(__name__)


def get(session: Session, *, id: int) -> Optional[Acquisition]:
    return session.query(Acquisition).filter(Acquisition.id == id).first()


def get_multi(session: Session, *, skip: int = 0, limit: int = 100) -> List[Optional[Acquisition]]:
    return session.query(Acquisition).offset(skip).limit(limit).all()


def create(session: Session, *, params: AcquisitionCreateModel) -> Acquisition:
    data = jsonable_encoder(params)
    entity = Acquisition(**data)
    session.add(entity)
    session.commit()
    session.refresh(entity)

    entity.location = os.path.join(
        entity.roi.acquisitions_location,
        ACQUISITION_LOCATION_FORMAT.format(id=entity.id),
    )
    if not os.path.exists(entity.location):
        logger.debug(f'Create location for acquisition {entity.id}: {entity.location}')
        os.makedirs(entity.location)

    session.commit()
    session.refresh(entity)
    return entity


def remove(session: Session, *, id: int):
    item = session.query(Acquisition).filter(Acquisition.id == id).first()
    session.delete(item)
    session.commit()
    return item
