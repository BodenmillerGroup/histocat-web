import logging
import os
from typing import List, Optional

from fastapi.encoders import jsonable_encoder
from sqlalchemy.orm import Session
from sqlalchemy.orm.attributes import flag_modified

from .models import DATASET_LOCATION_FORMAT, DatasetModel
from .dto import DatasetCreateDto, DatasetUpdateDto

logger = logging.getLogger(__name__)


def get(session: Session, *, id: int) -> Optional[DatasetModel]:
    return session.query(DatasetModel).filter(DatasetModel.id == id).first()


def get_own_by_experiment_id(session: Session, *, experiment_id: int) -> List[DatasetModel]:
    return session.query(DatasetModel).filter(DatasetModel.experiment_id == experiment_id).all()


def get_multi(session: Session, *, skip: int = 0, limit: int = 1000) -> List[DatasetModel]:
    return session.query(DatasetModel).offset(skip).limit(limit).all()


def create(session: Session, *, params: DatasetCreateDto) -> DatasetModel:
    data = jsonable_encoder(params)
    entity = DatasetModel(**data)
    session.add(entity)
    session.commit()
    session.refresh(entity)

    entity.location = os.path.join(entity.experiment.datasets_location, DATASET_LOCATION_FORMAT.format(id=entity.id),)
    if not os.path.exists(entity.location):
        logger.debug(f"Create location for dataset {entity.id}: {entity.location}")
        os.makedirs(entity.location)

    session.commit()
    session.refresh(entity)
    return entity


def update(session: Session, *, item: DatasetModel, params: DatasetUpdateDto) -> DatasetModel:
    data = jsonable_encoder(item)
    update_data = params.dict(skip_defaults=True)
    for field in data:
        if field in update_data:
            setattr(item, field, update_data[field])
    session.add(item)
    session.commit()
    session.refresh(item)
    return item


def update_output(session: Session, *, dataset_id: int, result_type: str, result: dict) -> DatasetModel:
    item = session.query(DatasetModel).filter(DatasetModel.id == dataset_id).first()
    session.refresh(item, attribute_names=["output"])

    output = item.output if item.output else {}
    result_output = output.get(result_type) if result_type in output else {}
    result_output[result.get("name")] = result
    output[result_type] = result_output
    item.output = output

    # TODO: https://stackoverflow.com/questions/42559434/updates-to-json-field-dont-persist-to-db
    flag_modified(item, "output")
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
