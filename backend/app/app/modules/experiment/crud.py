import logging
import os
from typing import List, Optional, Set

from fastapi.encoders import jsonable_encoder
from sqlalchemy import func
from sqlalchemy.orm import Session

from app.core.config import ROOT_DATA_DIRECTORY
from .db import Experiment, EXPERIMENT_LOCATION_FORMAT
from .models import ExperimentCreateModel, ExperimentUpdateModel

logger = logging.getLogger(__name__)


def get(session: Session, *, id: int) -> Optional[Experiment]:
    return session.query(Experiment).filter(Experiment.id == id).first()


def get_by_name(session: Session, *, name: str) -> Optional[Experiment]:
    return session.query(Experiment).filter(Experiment.name == name).first()


def get_multi(session: Session, *, skip: int = 0, limit: int = 100) -> List[Optional[Experiment]]:
    items = session.query(Experiment).offset(skip).limit(limit).all()
    return items


def get_tags(session: Session) -> Set[str]:
    items = session.query(func.unnest(Experiment.tags)).distinct().all()
    return {e[0] for e in items}


def create(session: Session, *, user_id: int, params: ExperimentCreateModel) -> Experiment:
    entity = Experiment(
        user_id=user_id,
        name=params.name,
        description=params.description,
        meta=params.meta,
        tags=params.tags
    )
    session.add(entity)
    session.commit()
    session.refresh(entity)

    entity.location = os.path.join(
        ROOT_DATA_DIRECTORY,
        EXPERIMENT_LOCATION_FORMAT.format(id=entity.id),
    )
    if not os.path.exists(entity.location):
        logger.debug(f'Create location for experiment {entity.id}: {entity.location}')
        os.makedirs(entity.location)

    session.commit()
    session.refresh(entity)

    return entity


def update(session: Session, *, item: Experiment, params: ExperimentUpdateModel) -> Experiment:
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
    item = session.query(Experiment).filter(Experiment.id == id).first()
    session.delete(item)
    session.commit()
    return item


def get_data(session: Session, *, id: int) -> Optional[Experiment]:
    return session.query(Experiment).filter(Experiment.id == id).first()
