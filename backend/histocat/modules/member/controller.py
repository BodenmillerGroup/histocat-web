from fastapi import APIRouter, Depends, HTTPException
from sqlalchemy.orm import Session

from histocat.api.db import get_db
from histocat.api.security import get_active_user, get_admin
from histocat.modules.user.models import UserModel

from . import service
from .dto import MemberCreateDto, MemberDto, MemberUpdateDto

router = APIRouter()


@router.post("/members", response_model=MemberDto)
def create_group_member(
    params: MemberCreateDto, db: Session = Depends(get_db), user: UserModel = Depends(get_active_user),
):
    """Create new member."""
    item = service.get_by_group_id_and_user_id(db, group_id=params.group_id, user_id=params.user_id)
    if item:
        raise HTTPException(
            status_code=400, detail="Member already exists.",
        )
    item = service.create(db, params=params)
    return item


@router.get("/members/{member_id}", response_model=MemberDto)
def get_by_id(member_id: int, db: Session = Depends(get_db), user: UserModel = Depends(get_active_user)):
    """Get member by id."""
    item = service.get_by_id(db, member_id)
    return item


@router.patch("/members/{member_id}", response_model=MemberDto)
def update(
    member_id: int, params: MemberUpdateDto, db: Session = Depends(get_db), user: UserModel = Depends(get_admin),
):
    """Update the member."""
    item = service.get_by_id(db, member_id)
    item = service.update(db, item=item, params=params)
    return item


@router.delete("/members/{member_id}", response_model=MemberDto)
def delete_by_id(
    member_id: int, user: UserModel = Depends(get_active_user), db: Session = Depends(get_db),
):
    """
    Delete member by id
    """
    item = service.delete_by_id(db, id=member_id)
    return item
