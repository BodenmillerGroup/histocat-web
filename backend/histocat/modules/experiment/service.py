import logging
import os
from typing import List, Optional, Set

from sqlalchemy import func
from sqlalchemy.orm import Session

from histocat.config import config

from .dto import ExperimentCreateDto, ExperimentUpdateDto
from .models import EXPERIMENT_LOCATION_FORMAT, ExperimentModel

logger = logging.getLogger(__name__)


def get(session: Session, *, id: int) -> Optional[ExperimentModel]:
    return session.query(ExperimentModel).filter(ExperimentModel.id == id).first()


def get_by_name(session: Session, *, name: str) -> Optional[ExperimentModel]:
    return session.query(ExperimentModel).filter(ExperimentModel.name == name).first()


def get_group_experiments(session: Session, *, group_id: int) -> List[Optional[ExperimentModel]]:
    return session.query(ExperimentModel).filter(ExperimentModel.group_id == group_id).all()


def get_tags(session: Session) -> Set[str]:
    items = session.query(func.unnest(ExperimentModel.tags)).distinct().all()
    return {e[0] for e in items}


def create(session: Session, *, group_id: int, params: ExperimentCreateDto) -> ExperimentModel:
    entity = ExperimentModel(group_id=group_id, name=params.name, description=params.description, tags=params.tags)
    session.add(entity)
    session.commit()
    session.refresh(entity)

    entity.location = os.path.join(config.ROOT_DATA_DIRECTORY, EXPERIMENT_LOCATION_FORMAT.format(id=entity.id))
    if not os.path.exists(entity.location):
        logger.debug(f"Create location for experiment {entity.id}: {entity.location}")
        os.makedirs(entity.location)

    session.commit()
    session.refresh(entity)

    return entity


def update(session: Session, *, item: ExperimentModel, params: ExperimentUpdateDto) -> ExperimentModel:
    data = item.as_dict()
    update_data = params.dict(exclude_unset=True)
    for field in data:
        if field in update_data:
            setattr(item, field, update_data[field])
    session.add(item)
    session.commit()
    session.refresh(item)
    return item


def remove(session: Session, *, id: int):
    item = session.query(ExperimentModel).filter(ExperimentModel.id == id).first()
    session.delete(item)
    session.commit()
    return item


def get_data(session: Session, *, id: int) -> Optional[ExperimentModel]:
    result = session.query(ExperimentModel).filter(ExperimentModel.id == id).first()
    return result
