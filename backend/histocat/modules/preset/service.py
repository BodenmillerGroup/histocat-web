import logging
from typing import Optional, Sequence

from fastapi.encoders import jsonable_encoder
from sqlalchemy.orm import Session

from .dto import PresetCreateDto
from .models import PresetModel

logger = logging.getLogger(__name__)


def get_by_id(session: Session, *, id: int) -> Optional[PresetModel]:
    return session.query(PresetModel).filter(PresetModel.id == id).first()


def get_project_presets(session: Session, *, project_id: int) -> Sequence[PresetModel]:
    return (
        session.query(PresetModel.id, PresetModel.name, PresetModel.description, PresetModel.created_at)
        .filter(PresetModel.project_id == project_id)
        .all()
    )


def create(session: Session, *, params: PresetCreateDto, member_id: int) -> PresetModel:
    data = jsonable_encoder(params)
    data["member_id"] = member_id
    entity = PresetModel(**data)
    session.add(entity)
    session.commit()
    session.refresh(entity)
    return entity


def delete_by_id(session: Session, *, id: int) -> int:
    item = session.query(PresetModel).filter(PresetModel.id == id).first()
    session.delete(item)
    session.commit()
    return id
