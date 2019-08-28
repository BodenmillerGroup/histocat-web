import logging
from typing import List, Optional

from fastapi.encoders import jsonable_encoder
from sqlalchemy.orm import Session

from .db import AcquisitionArtifact
from .models import AcquisitionArtifactCreateModel

logger = logging.getLogger(__name__)


def get_by_id(session: Session, *, id: int) -> Optional[AcquisitionArtifact]:
    return session.query(AcquisitionArtifact).filter(AcquisitionArtifact.id == id).first()


def get_multi(session: Session, *, skip: int = 0, limit: int = 100) -> List[Optional[AcquisitionArtifact]]:
    return session.query(AcquisitionArtifact).offset(skip).limit(limit).all()


def create(session: Session, *, params: AcquisitionArtifactCreateModel) -> AcquisitionArtifact:
    data = jsonable_encoder(params)
    entity = AcquisitionArtifact(**data)
    session.add(entity)
    session.commit()
    session.refresh(entity)
    return entity


def remove_by_id(session: Session, *, id: int):
    item = session.query(AcquisitionArtifact).filter(AcquisitionArtifact.id == id).first()
    session.delete(item)
    session.commit()
    return item
