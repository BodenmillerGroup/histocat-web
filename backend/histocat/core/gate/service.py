import logging
from typing import Optional, Sequence

from sqlalchemy.orm import Session

from .dto import GateCreateDto, GateUpdateDto
from .models import GateModel

logger = logging.getLogger(__name__)


def get_by_id(session: Session, *, id: int) -> Optional[GateModel]:
    return session.query(GateModel).filter(GateModel.id == id).first()


def get_dataset_gates(session: Session, *, dataset_id: int) -> Sequence[GateModel]:
    return session.query(GateModel).filter(GateModel.dataset_id == dataset_id).order_by(GateModel.id).all()


def create(session: Session, *, params: GateCreateDto) -> GateModel:
    data = params.dict()
    entity = GateModel(**data)
    session.add(entity)
    session.commit()
    session.refresh(entity)
    return entity


def update(session: Session, *, item: GateModel, params: GateUpdateDto) -> GateModel:
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
    item = session.query(GateModel).filter(GateModel.id == id).first()
    session.delete(item)
    session.commit()
    return id
