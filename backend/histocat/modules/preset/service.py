import logging
from typing import Sequence, Optional

from fastapi.encoders import jsonable_encoder
from sqlalchemy.orm import Session

from .dto import PresetCreateDto
from .models import PresetModel

logger = logging.getLogger(__name__)


def get_by_id(session: Session, *, id: int) -> Optional[PresetModel]:
    return session.query(PresetModel).filter(PresetModel.id == id).first()


def get_experiment_presets(session: Session, *, experiment_id: int) -> Sequence[PresetModel]:
    return session.query(PresetModel).filter(PresetModel.experiment_id == experiment_id).all()


def create(session: Session, *, params: PresetCreateDto) -> PresetModel:
    data = jsonable_encoder(params)
    entity = PresetModel(**data)
    session.add(entity)
    session.commit()
    session.refresh(entity)
    return entity


def delete_by_id(session: Session, *, id: int) -> Optional[PresetModel]:
    item = session.query(PresetModel).filter(PresetModel.id == id).first()
    session.delete(item)
    session.commit()
    return item
