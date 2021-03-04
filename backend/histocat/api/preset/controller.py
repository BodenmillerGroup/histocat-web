from typing import Sequence

from fastapi import APIRouter, Depends
from sqlalchemy.orm import Session

from histocat.api.db import get_db
from histocat.api.security import get_active_member
from histocat.core.member.models import MemberModel
from histocat.core.preset import service
from histocat.core.preset.dto import PresetCreateDto, PresetDto

router = APIRouter()


@router.get("/groups/{group_id}/presets/{preset_id}", response_model=PresetDto)
def get_by_id(
    group_id: int, preset_id: int, member: MemberModel = Depends(get_active_member), db: Session = Depends(get_db),
):
    """
    Get preset by id
    """
    item = service.get_by_id(db, id=preset_id)
    return item


@router.get("/groups/{group_id}/projects/{project_id}/presets", response_model=Sequence[PresetDto])
def get_project_presets(
    group_id: int, project_id: int, member: MemberModel = Depends(get_active_member), db: Session = Depends(get_db),
):
    """
    Get all project presets
    """
    item = service.get_project_presets(db, project_id=project_id)
    return item


@router.post("/groups/{group_id}/presets", response_model=PresetDto)
def create(
    group_id: int,
    params: PresetCreateDto,
    db: Session = Depends(get_db),
    member: MemberModel = Depends(get_active_member),
):
    """
    Create new preset
    """
    items = service.create(db, params=params, member_id=member.id)
    return items


@router.delete("/groups/{group_id}/presets/{preset_id}", response_model=int)
def delete_by_id(
    group_id: int, preset_id: int, member: MemberModel = Depends(get_active_member), db: Session = Depends(get_db),
):
    """
    Delete preset by id
    """
    service.delete_by_id(db, id=preset_id)
    return preset_id
