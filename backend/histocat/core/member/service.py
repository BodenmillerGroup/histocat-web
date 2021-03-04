from typing import Optional, Sequence

from sqlalchemy.orm import Session

from .dto import MemberCreateDto, MemberUpdateDto
from .models import MemberModel


def get_by_id(session: Session, id: int) -> Optional[MemberModel]:
    return session.query(MemberModel).filter(MemberModel.id == id).first()


def get_by_group_id_and_user_id(session: Session, group_id: int, user_id: int) -> Optional[MemberModel]:
    return session.query(MemberModel).filter(MemberModel.group_id == group_id, MemberModel.user_id == user_id).first()


def get_group_members(session: Session, *, group_id: int) -> Sequence[MemberModel]:
    return session.query(MemberModel).filter(MemberModel.group_id == group_id).all()


def create(session: Session, *, group_id: int, params: MemberCreateDto) -> MemberModel:
    entity = MemberModel(group_id=group_id, user_id=params.user_id, role=params.role, is_active=params.is_active,)
    session.add(entity)
    session.commit()
    session.refresh(entity)
    return entity


def update(session: Session, *, item: MemberModel, params: MemberUpdateDto) -> MemberModel:
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
    item = session.query(MemberModel).filter(MemberModel.id == id).first()
    session.delete(item)
    session.commit()
    return id
