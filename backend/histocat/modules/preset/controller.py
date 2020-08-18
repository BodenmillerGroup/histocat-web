from typing import Sequence

from fastapi import APIRouter, Depends
from sqlalchemy.orm import Session

from histocat.api.utils.db import get_db
from histocat.api.utils.security import get_current_active_user
from histocat.modules.user.models import UserModel

from . import service
from .dto import PresetCreateDto, PresetDto

router = APIRouter()


@router.get("/presets/{preset_id}", response_model=PresetDto)
def get_by_id(
    preset_id: int, current_user: UserModel = Depends(get_current_active_user), db: Session = Depends(get_db),
):
    """
    Get preset by id
    """
    item = service.get_by_id(db, id=preset_id)
    return item


@router.get("/experiments/{experiment_id}/presets", response_model=Sequence[PresetDto])
def get_experiment_presets(
    experiment_id: int, current_user: UserModel = Depends(get_current_active_user), db: Session = Depends(get_db),
):
    """
    Get all experiment presets
    """
    item = service.get_experiment_presets(db, experiment_id=experiment_id)
    return item


@router.post("/presets", response_model=PresetDto)
def create(
    *,
    db: Session = Depends(get_db),
    params: PresetCreateDto,
    current_user: UserModel = Depends(get_current_active_user),
):
    """
    Create new preset
    """
    items = service.create(db, params=params)
    return items


@router.delete("/presets/{preset_id}", response_model=int)
def delete_by_id(
    preset_id: int, current_user: UserModel = Depends(get_current_active_user), db: Session = Depends(get_db),
):
    """
    Delete preset by id
    """
    service.delete_by_id(db, id=preset_id)
    return preset_id
