import logging
import os
from typing import Optional, Sequence

from sqlalchemy.orm import Session

from histocat.config import config

from .dto import ModelCreateDto, ModelUpdateDto
from .models import ModelModel

logger = logging.getLogger(__name__)


def get(session: Session, *, id: int) -> Optional[ModelModel]:
    return session.query(ModelModel).filter(ModelModel.id == id).first()


def get_by_name(session: Session, *, name: str) -> Optional[ModelModel]:
    return session.query(ModelModel).filter(ModelModel.name == name).first()


def get_all_models(session: Session) -> Sequence[Optional[ModelModel]]:
    return session.query(ModelModel).order_by(ModelModel.id.desc()).all()


def create(session: Session, *, params: ModelCreateDto) -> ModelModel:
    entity = ModelModel(name=params.name, description=params.description, application=params.application)
    session.add(entity)
    session.commit()
    session.refresh(entity)

    entity.location = os.path.join(config.MODELS_DIRECTORY, str(entity.id))
    if not os.path.exists(entity.location):
        logger.debug(f"Create location for model {entity.id}: {entity.location}")
        os.makedirs(entity.location)

    session.commit()
    session.refresh(entity)

    return entity


def update(session: Session, *, item: ModelModel, params: ModelUpdateDto) -> ModelModel:
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
    item = session.query(ModelModel).filter(ModelModel.id == id).first()
    session.delete(item)
    session.commit()
    return item
