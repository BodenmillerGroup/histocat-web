from typing import Sequence, Set

from fastapi import APIRouter, Depends, HTTPException
from sqlalchemy.orm import Session

from histocat.api.utils.db import get_db
from histocat.api.utils.security import (
    get_current_active_superuser,
    get_current_active_user,
)
from histocat.modules.user.models import UserModel

from . import service
from .dto import GroupCreateDto, GroupDto, GroupUpdateDto

router = APIRouter()


@router.get("/groups/tags", response_model=Set[str])
def get_tags(db: Session = Depends(get_db), current_user: UserModel = Depends(get_current_active_user)):
    """
    Get groups tags
    """
    return service.get_tags(db)


@router.get("/groups", response_model=Sequence[GroupDto])
def get_all(db: Session = Depends(get_db), current_user: UserModel = Depends(get_current_active_user)):
    """Get all groups."""
    items = service.get_all(db, user_id=current_user.id)
    return items


@router.post("/groups", response_model=GroupDto)
def create(
    params: GroupCreateDto,
    db: Session = Depends(get_db),
    current_user: UserModel = Depends(get_current_active_user),
):
    """Create new group."""
    item = service.get_by_name(db, name=params.name)
    if item:
        raise HTTPException(
            status_code=400, detail="The group with this name already exists.",
        )
    item = service.create(db, params=params, user_id=current_user.id)
    return item


@router.get("/groups/{id}", response_model=GroupDto)
def get_by_id(id: int, db: Session = Depends(get_db), current_user: UserModel = Depends(get_current_active_user)):
    """Get group by id."""
    item = service.get_by_id(db, id)
    return item


@router.patch("/groups/{id}", response_model=GroupDto)
def update(
    id: int,
    params: GroupUpdateDto,
    db: Session = Depends(get_db),
    current_user: UserModel = Depends(get_current_active_superuser),
):
    """Update the group."""
    item = service.get_by_id(db, id)
    item = service.update(db, item=item, params=params)
    return item


@router.post("/groups/{id}/join", response_model=GroupDto)
def join(
    id: int,
    db: Session = Depends(get_db),
    current_user: UserModel = Depends(get_current_active_user),
):
    """Join the group."""
    return service.join(db, group_id=id, user_id=current_user.id)


@router.delete("/groups/{id}", response_model=int)
def delete_by_id(
    id: int, current_user: UserModel = Depends(get_current_active_user), db: Session = Depends(get_db),
):
    """
    Delete group by id
    """
    item = service.delete_by_id(db, id=id)
    return item.id
