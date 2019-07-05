import logging
import os
from typing import List, Optional

from sqlalchemy.orm import Session

from .db import Dataset, DATASET_LOCATION_FORMAT
from .models import DatasetCreateModel

logger = logging.getLogger(__name__)


def get(session: Session, *, id: int) -> Optional[Dataset]:
    return session.query(Dataset).filter(Dataset.id == id).first()


def get_own_by_experiment_id(session: Session, *, user_id: int, experiment_id: int) -> List[Dataset]:
    return session.query(Dataset).filter(Dataset.user_id == user_id, Dataset.experiment_id == experiment_id).all()


def get_multi(session: Session, *, skip: int = 0, limit: int = 100) -> List[Dataset]:
    return session.query(Dataset).offset(skip).limit(limit).all()


def create(session: Session, *, user_id: int, params: DatasetCreateModel) -> Dataset:
    entity = Dataset(
        user_id=user_id,
        experiment_id=params.experiment_id,
        name=params.name,
        description=params.description,
        meta=params.meta,
    )
    session.add(entity)
    session.commit()
    session.refresh(entity)

    entity.location = os.path.join(
        entity.experiment.datasets_location,
        DATASET_LOCATION_FORMAT.format(id=entity.id),
    )
    if not os.path.exists(entity.location):
        logger.debug(f'Create location for dataset {entity.id}: {entity.location}')
        os.makedirs(entity.location)

    session.commit()
    session.refresh(entity)
    return entity


def remove(session: Session, *, id: int):
    item = session.query(Dataset).filter(Dataset.id == id).first()
    if item:
        session.delete(item)
        session.commit()
        return item
