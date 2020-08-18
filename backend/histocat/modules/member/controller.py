from fastapi import APIRouter, Depends, HTTPException
from sqlalchemy.orm import Session

from histocat.api.utils.db import get_db
from histocat.api.utils.security import (
    get_current_active_superuser,
    get_current_active_user,
)
from histocat.modules.user.models import UserModel

from . import service
from .dto import MemberCreateDto, MemberDto, MemberUpdateDto

router = APIRouter()


@router.post("/", response_model=MemberDto)
def create(
    params: MemberCreateDto,
    db: Session = Depends(get_db),
    current_user: UserModel = Depends(get_current_active_superuser),
):
    """Create new member."""
    item = service.get_by_group_id_and_user_id(db, group_id=params.group_id, user_id=params.user_id)
    if item:
        raise HTTPException(
            status_code=400, detail="Member already exists.",
        )
    item = service.create(db, params=params)
    return item


@router.get("/{id}", response_model=MemberDto)
def get_by_id(id: int, db: Session = Depends(get_db), current_user: UserModel = Depends(get_current_active_user)):
    """Get member by id."""
    item = service.get_by_id(db, id)
    return item


@router.put("/{id}", response_model=MemberDto)
def update(
    id: int,
    params: MemberUpdateDto,
    db: Session = Depends(get_db),
    current_user: UserModel = Depends(get_current_active_superuser),
):
    """Update the member."""
    item = service.get_by_id(db, id)
    item = service.update(db, item=item, params=params)
    return item


@router.delete("/{id}", response_model=MemberDto)
def delete_by_id(
    id: int, current_user: UserModel = Depends(get_current_active_user), db: Session = Depends(get_db),
):
    """
    Delete member by id
    """
    item = service.delete_by_id(db, id=id)
    return item
