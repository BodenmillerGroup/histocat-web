from typing import Sequence

from fastapi import APIRouter, Depends
from sqlalchemy.orm import Session

from histocat.api.db import get_db
from histocat.api.security import get_active_user
from histocat.modules.user.models import UserModel

from . import service
from .dto import PresetCreateDto, PresetDto

router = APIRouter()


@router.get("/presets/{preset_id}", response_model=PresetDto)
def get_by_id(
    preset_id: int, user: UserModel = Depends(get_active_user), db: Session = Depends(get_db),
):
    """
    Get preset by id
    """
    item = service.get_by_id(db, id=preset_id)
    return item


@router.get("/projects/{project_id}/presets", response_model=Sequence[PresetDto])
def get_project_presets(
    project_id: int, user: UserModel = Depends(get_active_user), db: Session = Depends(get_db),
):
    """
    Get all project presets
    """
    item = service.get_project_presets(db, project_id=project_id)
    return item


@router.post("/presets", response_model=PresetDto)
def create(
    *, db: Session = Depends(get_db), params: PresetCreateDto, user: UserModel = Depends(get_active_user),
):
    """
    Create new preset
    """
    items = service.create(db, params=params)
    return items


@router.delete("/presets/{preset_id}", response_model=int)
def delete_by_id(
    preset_id: int, user: UserModel = Depends(get_active_user), db: Session = Depends(get_db),
):
    """
    Delete preset by id
    """
    service.delete_by_id(db, id=preset_id)
    return preset_id
