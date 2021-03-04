import logging
import os
from typing import Optional, Sequence, Set

from sqlalchemy import func
from sqlalchemy.orm import Session

from histocat.config import config
from histocat.core.member import service as member_service
from histocat.core.member.dto import MemberCreateDto
from histocat.core.member.models import MemberModel

from .dto import GroupCreateDto, GroupUpdateDto
from .models import GroupModel

logger = logging.getLogger(__name__)


def get_by_id(session: Session, id: int) -> Optional[GroupModel]:
    return session.query(GroupModel).filter(GroupModel.id == id).first()


def get_by_name(session: Session, *, name: str) -> Optional[GroupModel]:
    return session.query(GroupModel).filter(GroupModel.name == name).first()


def get_all(session: Session, *, user_id: int) -> Sequence[GroupModel]:
    return (
        session.query(GroupModel)
        .outerjoin(GroupModel.members)
        .filter(GroupModel.is_open | (MemberModel.user_id == user_id))
        .all()
    )


def get_tags(session: Session) -> Set[str]:
    items = session.query(func.unnest(GroupModel.tags)).distinct().all()
    return {e[0] for e in items}


def create(session: Session, *, params: GroupCreateDto, user_id: int) -> GroupModel:
    entity = GroupModel(
        name=params.name, description=params.description, url=params.url, is_open=params.is_open, tags=params.tags
    )
    session.add(entity)
    session.commit()
    session.refresh(entity)

    entity.location = os.path.join(config.ROOT_DATA_DIRECTORY, str(entity.id))
    if not os.path.exists(entity.location):
        logger.debug(f"Create location for group {entity.id}: {entity.location}")
        os.makedirs(entity.location)

    session.commit()
    session.refresh(entity)

    member = member_service.create(
        session, group_id=entity.id, params=MemberCreateDto(user_id=user_id, role=100, is_active=True)
    )

    return entity


def update(session: Session, *, item: GroupModel, params: GroupUpdateDto) -> GroupModel:
    data = item.as_dict()
    update_data = params.dict(exclude_unset=True)
    for field in data:
        if field in update_data:
            setattr(item, field, update_data[field])
    session.add(item)
    session.commit()
    session.refresh(item)
    return item


def join(session: Session, *, group_id: int, user_id: int) -> Optional[GroupModel]:
    member = member_service.create(
        session, group_id=group_id, params=MemberCreateDto(user_id=user_id, role=10, is_active=True)
    )
    return get_by_id(session, id=group_id)


def delete_by_id(session: Session, *, id: int):
    item = session.query(GroupModel).filter(GroupModel.id == id).first()
    session.delete(item)
    session.commit()
    return item
