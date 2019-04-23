from typing import List

from fastapi import APIRouter, Depends, HTTPException
from sqlalchemy.orm import Session

from app.api.utils.db import get_db
from app.api.utils.security import get_current_active_superuser, get_current_active_user
from app.modules.user.db import User
from . import crud
from .models import ChannelModel, ChannelCreateModel, ChannelUpdateModel

router = APIRouter()


@router.get("/", response_model=List[ChannelModel])
def read_channels(
    db: Session = Depends(get_db),
    skip: int = 0,
    limit: int = 100,
    current_user: User = Depends(get_current_active_superuser),
):
    """
    Retrieve channels
    """
    items = crud.get_multi(db, skip=skip, limit=limit)
    return items


@router.post("/", response_model=ChannelModel)
def create_channel(
    *,
    db: Session = Depends(get_db),
    params: ChannelCreateModel,
    current_user: User = Depends(get_current_active_superuser),
):
    """
    Create new channel
    """
    item = crud.get_by_name(db, name=params.name)
    if item:
        raise HTTPException(
            status_code=400,
            detail="The channel with this name already exists in the system.",
        )
    item = crud.create(db, params=params)
    return item


@router.get("/{id}", response_model=ChannelModel)
def read_channel_by_id(
    id: int,
    current_user: User = Depends(get_current_active_user),
    db: Session = Depends(get_db),
):
    """
    Get a specific channel by id
    """
    item = crud.get(db, id=id)
    return item


@router.put("/{id}", response_model=ChannelModel)
def update_channel(
    *,
    db: Session = Depends(get_db),
    id: int,
    params: ChannelUpdateModel,
    current_user: User = Depends(get_current_active_superuser),
):
    """
    Update a channel
    """
    item = crud.get(db, id=id)

    if not item:
        raise HTTPException(
            status_code=404,
            detail="The channel with this id does not exist in the system",
        )
    item = crud.update(db, item=item, params=params)
    return item
