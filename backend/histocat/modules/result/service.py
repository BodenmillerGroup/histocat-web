import logging
import os
from typing import List, Optional

from fastapi.encoders import jsonable_encoder
from sqlalchemy.orm import Session

from .dto import ResultCreateDto, ResultUpdateDto
from .models import ResultModel

logger = logging.getLogger(__name__)


def get(session: Session, *, id: int) -> Optional[ResultModel]:
    return session.query(ResultModel).filter(ResultModel.id == id).first()


def get_dataset_results(session: Session, *, dataset_id: int) -> List[ResultModel]:
    return session.query(ResultModel).filter(ResultModel.dataset_id == dataset_id, ResultModel.status == "ready").all()


def create(session: Session, *, params: ResultCreateDto) -> ResultModel:
    data = jsonable_encoder(params)
    entity = ResultModel(**data)
    session.add(entity)
    session.commit()
    session.refresh(entity)

    entity.location = os.path.join(entity.dataset.results_location, str(entity.id))
    if not os.path.exists(entity.location):
        logger.debug(f"Create location for result {entity.id}: {entity.location}")
        os.makedirs(entity.location)

    session.commit()
    session.refresh(entity)
    return entity


def update(session: Session, *, item: ResultModel, params: ResultUpdateDto) -> ResultModel:
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
    item = session.query(ResultModel).filter(ResultModel.id == id).first()
    if item:
        session.delete(item)
        session.commit()
        return item
