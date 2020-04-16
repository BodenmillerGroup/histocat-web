from typing import List

from fastapi import APIRouter, Depends
from fastapi.responses import ORJSONResponse
from sqlalchemy.orm import Session

from app.api.utils.db import get_db
from app.api.utils.security import get_current_active_user
from app.modules.user.models import UserModel

from . import service
from .dto import ShareCreateDto, ShareDto

router = APIRouter()


@router.get("/{experiment_id}", response_model=List[ShareDto])
def read_all_by_experiment_id(
    experiment_id: int, current_user: UserModel = Depends(get_current_active_user), db: Session = Depends(get_db),
):
    """
    Retrieve all experiment shares
    """
    items = service.get_by_experiment_id(db, experiment_id=experiment_id)
    return ORJSONResponse(items)


@router.post("/", response_model=List[ShareDto])
def create(
    *,
    db: Session = Depends(get_db),
    params: ShareCreateDto,
    current_user: UserModel = Depends(get_current_active_user),
):
    """
    Create new share
    """
    items = service.create(db, params=params)
    return ORJSONResponse(items)


@router.delete("/{user_id}/{experiment_id}", response_model=ShareDto)
def delete(
    user_id: int,
    experiment_id: int,
    current_user: UserModel = Depends(get_current_active_user),
    db: Session = Depends(get_db),
):
    """
    Delete share
    """
    item = service.remove(db, user_id=user_id, experiment_id=experiment_id)
    return ORJSONResponse(item)
