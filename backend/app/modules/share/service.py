import logging
from typing import List, Optional

from sqlalchemy.orm import Session

from .dto import ShareCreateDto
from .models import ShareModel

logger = logging.getLogger(__name__)


def get_by_user_id(session: Session, *, user_id: int) -> Optional[List[ShareModel]]:
    return session.query(ShareModel).filter(ShareModel.user_id == user_id).all()


def get_by_experiment_id(session: Session, *, experiment_id: int) -> Optional[List[ShareModel]]:
    return session.query(ShareModel).filter(ShareModel.experiment_id == experiment_id).all()


def create(session: Session, *, params: ShareCreateDto) -> List[ShareModel]:
    shares = []
    items = session.query(ShareModel).filter(ShareModel.experiment_id == params.experiment_id).all()
    for item in items:
        session.delete(item)
        session.commit()
    for user_id in params.user_ids:
        entity = ShareModel(user_id=user_id, experiment_id=params.experiment_id, permissions=params.permissions,)
        session.add(entity)
        session.commit()
        session.refresh(entity)
        shares.append(entity)
    return shares


def remove(session: Session, *, user_id: int, experiment_id: int):
    item = (
        session.query(ShareModel)
        .filter(ShareModel.user_id == user_id, ShareModel.experiment_id == experiment_id)
        .first()
    )
    session.delete(item)
    session.commit()
    return item
