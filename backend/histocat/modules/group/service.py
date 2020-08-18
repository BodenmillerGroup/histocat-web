from typing import Optional, Sequence

from sqlalchemy.orm import Session
from sqlalchemy import or_, and_, any_

from .dto import GroupCreateDto, GroupUpdateDto
from .models import GroupModel
from histocat.modules.member import service as member_service
from histocat.modules.member.dto import MemberCreateDto
from histocat.modules.member.models import MemberModel


def get_by_id(session: Session, id: int) -> Optional[GroupModel]:
    return session.query(GroupModel).filter(GroupModel.id == id).first()


def get_by_name(session: Session, *, name: str) -> Optional[GroupModel]:
    return session.query(GroupModel).filter(GroupModel.name == name).first()


def get_all(session: Session, *, user_id: int) -> Sequence[GroupModel]:
    return session.query(GroupModel).outerjoin(GroupModel.members).filter((GroupModel.is_open) | (MemberModel.user_id == user_id)).all()


def create(session: Session, *, params: GroupCreateDto, user_id: int) -> GroupModel:
    entity = GroupModel(name=params.name, description=params.description, url=params.url, is_open=params.is_open)
    session.add(entity)
    session.commit()
    session.refresh(entity)

    member = member_service.create(session, params=MemberCreateDto(group_id=entity.id, user_id=user_id, role=100, is_active=True))

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


def delete_by_id(session: Session, *, id: int):
    item = session.query(GroupModel).filter(GroupModel.id == id).first()
    session.delete(item)
    session.commit()
    return item
