from fastapi import APIRouter, Depends
from sqlalchemy.orm import Session
from typing import List

from app.api.utils.db import get_db
from app.api.utils.security import get_current_active_superuser, get_current_active_user
from app.modules.user.db import User
from . import crud
from .models import ShareModel, ShareCreateModel

router = APIRouter()


@router.get("/{experiment_id}", response_model=List[ShareModel])
def read_all_by_experiment_id(
    experiment_id: int,
    current_user: User = Depends(get_current_active_user),
    db: Session = Depends(get_db),
):
    """
    Retrieve all experiment shares
    """
    items = crud.get_by_experiment_id(db, experiment_id=experiment_id)
    return items


@router.post("/", response_model=List[ShareModel])
def create(
    *,
    db: Session = Depends(get_db),
    params: ShareCreateModel,
    current_user: User = Depends(get_current_active_user),
):
    """
    Create new share
    """
    items = crud.create(db, params=params)
    return items


@router.delete("/{user_id}/{experiment_id}", response_model=ShareModel)
def delete(
    user_id: int,
    experiment_id: int,
    current_user: User = Depends(get_current_active_user),
    db: Session = Depends(get_db),
):
    """
    Delete share
    """
    item = crud.remove(db, user_id=user_id, experiment_id=experiment_id)
    return item
