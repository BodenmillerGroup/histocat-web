from fastapi import APIRouter, Depends, HTTPException
from sqlalchemy.orm import Session

from histocat.api.db import get_db
from histocat.api.security import get_active_member
from histocat.core.member import service
from histocat.core.member.dto import MemberCreateDto, MemberDto, MemberUpdateDto
from histocat.core.member.models import MemberModel

router = APIRouter()


@router.post("/groups/{group_id}/members", response_model=MemberDto)
def create_group_member(
    group_id: int,
    params: MemberCreateDto,
    db: Session = Depends(get_db),
    member: MemberModel = Depends(get_active_member),
):
    """Create new member."""
    item = service.get_by_group_id_and_user_id(db, group_id=group_id, user_id=params.user_id)
    if item:
        raise HTTPException(
            status_code=400, detail="Member already exists.",
        )
    item = service.create(db, group_id=group_id, params=params)
    return item


@router.get("/groups/{group_id}/members/{member_id}", response_model=MemberDto)
def get_by_id(
    group_id: int, member_id: int, db: Session = Depends(get_db), member: MemberModel = Depends(get_active_member)
):
    """Get member by id."""
    item = service.get_by_id(db, member_id)
    return item


@router.patch("/groups/{group_id}/members/{member_id}", response_model=MemberDto)
def update(
    group_id: int,
    member_id: int,
    params: MemberUpdateDto,
    db: Session = Depends(get_db),
    member: MemberModel = Depends(get_active_member),
):
    """Update the member."""
    item = service.get_by_id(db, member_id)
    item = service.update(db, item=item, params=params)
    return item


@router.delete("/groups/{group_id}/members/{member_id}", response_model=int)
def delete_by_id(
    group_id: int, member_id: int, member: MemberModel = Depends(get_active_member), db: Session = Depends(get_db),
):
    """
    Delete member by id
    """
    item = service.delete_by_id(db, id=member_id)
    return item
