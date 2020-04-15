import logging
from typing import Optional

from fastapi.encoders import jsonable_encoder
from sqlalchemy.orm import Session

from .dto import AcquisitionCreateDto
from .models import AcquisitionModel

logger = logging.getLogger(__name__)


def get_by_id(session: Session, id: int) -> Optional[AcquisitionModel]:
    return session.query(AcquisitionModel).filter(AcquisitionModel.id == id).first()


def get_by_origin_id(session: Session, *, slide_id: int, origin_id: int) -> Optional[AcquisitionModel]:
    return (
        session.query(AcquisitionModel)
        .filter(AcquisitionModel.slide_id == slide_id, AcquisitionModel.origin_id == origin_id)
        .first()
    )


def create(session: Session, params: AcquisitionCreateDto) -> AcquisitionModel:
    data = jsonable_encoder(params)
    item = AcquisitionModel(**data)
    session.add(item)
    session.commit()
    session.refresh(item)
    return item


def delete_by_id(session: Session, id: int):
    item = get_by_id(session, id)
    session.delete(item)
    session.commit()
    return item
