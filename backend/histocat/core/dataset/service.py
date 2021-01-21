import logging
import os
from typing import List, Optional

from sqlalchemy.orm import Session

from .dto import DatasetCreateDto, DatasetUpdateDto
from .models import DatasetModel

logger = logging.getLogger(__name__)


def get(session: Session, *, id: int) -> Optional[DatasetModel]:
    return session.query(DatasetModel).filter(DatasetModel.id == id).first()


def get_project_datasets(session: Session, *, project_id: int) -> List[DatasetModel]:
    return session.query(DatasetModel).filter(DatasetModel.project_id == project_id).all()


def get_multi(session: Session, *, skip: int = 0, limit: int = 1000) -> List[DatasetModel]:
    return session.query(DatasetModel).offset(skip).limit(limit).all()


def create(session: Session, *, params: DatasetCreateDto) -> DatasetModel:
    data = params.dict()
    entity = DatasetModel(**data)
    session.add(entity)
    session.commit()
    session.refresh(entity)

    entity.location = os.path.join(entity.project.datasets_location, str(entity.id))
    if not os.path.exists(entity.location):
        logger.debug(f"Create location for dataset {entity.id}: {entity.location}")
        os.makedirs(entity.location)

    session.commit()
    session.refresh(entity)
    return entity


def update(session: Session, *, item: DatasetModel, params: DatasetUpdateDto) -> DatasetModel:
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
    item = session.query(DatasetModel).filter(DatasetModel.id == id).first()
    if item:
        session.delete(item)
        session.commit()
        return item
