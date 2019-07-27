import logging
from typing import Optional, List

from fastapi.encoders import jsonable_encoder
from sqlalchemy.orm import Session

from .db import Share
from .models import ShareCreateModel

logger = logging.getLogger(__name__)


def get_by_user_id(session: Session, *, user_id: int) -> Optional[List[Share]]:
    return session.query(Share).filter(Share.user_id == user_id).all()


def create(session: Session, *, params: ShareCreateModel) -> List[Share]:
    shares = []
    for user_id in params.user_ids:
        entity = Share(
            user_id=user_id,
            experiment_id=params.experiment_id,
            permissions=params.permissions
        )
        session.add(entity)
        session.commit()
        session.refresh(entity)
        shares.append(entity)
    return shares


def remove(session: Session, *, user_id: int, experiment_id: int):
    item = session.query(Share).filter(Share.user_id == user_id, Share.experiment_id == experiment_id).first()
    session.delete(item)
    session.commit()
    return item
