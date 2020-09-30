import logging
from typing import Optional, Sequence

from fastapi.encoders import jsonable_encoder
from sqlalchemy.orm import Session

from .dto import PipelineCreateDto, PipelineUpdateDto
from .models import PipelineModel

logger = logging.getLogger(__name__)


def get_by_id(session: Session, *, id: int) -> Optional[PipelineModel]:
    return session.query(PipelineModel).filter(PipelineModel.id == id).first()


def get_project_pipelines(session: Session, *, project_id: int) -> Sequence[PipelineModel]:
    return (
        session.query(PipelineModel.id, PipelineModel.name, PipelineModel.description, PipelineModel.created_at)
        .filter(PipelineModel.project_id == project_id)
        .all()
    )


def create(session: Session, *, params: PipelineCreateDto) -> PipelineModel:
    data = jsonable_encoder(params)
    entity = PipelineModel(**data)
    session.add(entity)
    session.commit()
    session.refresh(entity)
    return entity


def update(session: Session, *, item: PipelineModel, params: PipelineUpdateDto) -> PipelineModel:
    data = item.as_dict()
    update_data = params.dict(exclude_unset=True)
    for field in data:
        if field in update_data:
            setattr(item, field, update_data[field])
    session.add(item)
    session.commit()
    session.refresh(item)
    return item


def delete_by_id(session: Session, *, id: int) -> int:
    item = session.query(PipelineModel).filter(PipelineModel.id == id).first()
    session.delete(item)
    session.commit()
    return id
