import logging
from typing import Sequence, Optional

from fastapi.encoders import jsonable_encoder
from sqlalchemy.orm import Session

from .dto import GateCreateDto
from .models import GateModel

logger = logging.getLogger(__name__)


def get_by_id(session: Session, *, id: int) -> Optional[GateModel]:
    return session.query(GateModel).filter(GateModel.id == id).first()


def get_dataset_gates(session: Session, *, dataset_id: int) -> Sequence[GateModel]:
    return session.query(GateModel).filter(GateModel.dataset_id == dataset_id).all()


def create(session: Session, *, params: GateCreateDto) -> GateModel:
    data = jsonable_encoder(params)
    entity = GateModel(**data)
    session.add(entity)
    session.commit()
    session.refresh(entity)
    return entity


def delete_by_id(session: Session, *, id: int) -> Optional[GateModel]:
    item = session.query(GateModel).filter(GateModel.id == id).first()
    session.delete(item)
    session.commit()
    return item
