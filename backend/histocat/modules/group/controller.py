from typing import Sequence, Set

from fastapi import APIRouter, Depends, HTTPException
from sqlalchemy.orm import Session

from histocat.api.db import get_db
from histocat.api.security import get_active_member, get_active_user, get_group_admin
from histocat.modules.member.dto import MemberDto
from histocat.modules.member.models import MemberModel
from histocat.modules.member.service import (
    get_by_group_id_and_user_id,
    get_group_members,
)
from histocat.modules.user.models import UserModel

from . import service
from .dto import GroupCreateDto, GroupDto, GroupUpdateDto

router = APIRouter()


@router.get("/groups/tags", response_model=Set[str])
def get_tags(db: Session = Depends(get_db), user: UserModel = Depends(get_active_user)):
    """
    Get groups tags
    """
    return service.get_tags(db)


@router.get("/groups", response_model=Sequence[GroupDto])
def get_all(db: Session = Depends(get_db), user: UserModel = Depends(get_active_user)):
    """Get all groups."""
    items = service.get_all(db, user_id=user.id)
    return items


@router.post("/groups", response_model=GroupDto)
def create(
    params: GroupCreateDto, db: Session = Depends(get_db), user: UserModel = Depends(get_active_user),
):
    """Create new group."""
    item = service.get_by_name(db, name=params.name)
    if item:
        raise HTTPException(
            status_code=400, detail="The group with this name already exists.",
        )
    item = service.create(db, params=params, user_id=user.id)
    return item


@router.get("/groups/{group_id}", response_model=GroupDto)
def get_by_id(group_id: int, db: Session = Depends(get_db), member: MemberModel = Depends(get_active_member)):
    """Get group by id."""
    item = service.get_by_id(db, group_id)
    return item


@router.patch("/groups/{group_id}", response_model=GroupDto)
def update(
    group_id: int,
    params: GroupUpdateDto,
    db: Session = Depends(get_db),
    member: MemberModel = Depends(get_group_admin),
):
    """Update the group."""
    item = service.get_by_id(db, group_id)
    item = service.update(db, item=item, params=params)
    return item


@router.post("/groups/{group_id}/join", response_model=GroupDto)
def join(
    group_id: int, db: Session = Depends(get_db), user: UserModel = Depends(get_active_user),
):
    """Join the group."""
    return service.join(db, group_id=group_id, user_id=user.id)


@router.delete("/groups/{group_id}", response_model=int)
def delete_by_id(
    group_id: int, member: MemberModel = Depends(get_group_admin), db: Session = Depends(get_db),
):
    """
    Delete group by id
    """
    item = service.delete_by_id(db, id=group_id)
    return item.id


@router.get("/groups/{group_id}/members/me", response_model=MemberDto)
def get_own_member(group_id: int, db: Session = Depends(get_db), member: MemberModel = Depends(get_active_member)):
    """Get own member of the group."""
    item = get_by_group_id_and_user_id(db, group_id=group_id, user_id=member.user.id)
    return item


@router.get("/groups/{group_id}/members", response_model=Sequence[MemberDto])
def find_group_members(group_id: int, db: Session = Depends(get_db), member: MemberModel = Depends(get_active_member)):
    """Get group's members."""
    items = get_group_members(db, group_id=group_id)
    return items
