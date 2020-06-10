from typing import Sequence

from fastapi import APIRouter, Depends
from sqlalchemy.orm import Session

from histocat.api.utils.db import get_db
from histocat.api.utils.security import get_current_active_user
from histocat.modules.user.models import UserModel

from . import service
from .dto import GateDto, GateCreateDto

router = APIRouter()


@router.get("/gates/{gate_id}", response_model=GateDto)
def get_by_id(
    gate_id: int, current_user: UserModel = Depends(get_current_active_user), db: Session = Depends(get_db),
):
    """
    Get gate by id
    """
    item = service.get_by_id(db, id=gate_id)
    return item


@router.get("/datasets/{dataset_id}/gates", response_model=Sequence[GateDto])
def get_dataset_gates(
    dataset_id: int, current_user: UserModel = Depends(get_current_active_user), db: Session = Depends(get_db),
):
    """
    Get all dataset gates
    """
    item = service.get_dataset_gates(db, dataset_id=dataset_id)
    return item


@router.post("/gates", response_model=GateDto)
def create(
    *,
    db: Session = Depends(get_db),
    params: GateCreateDto,
    current_user: UserModel = Depends(get_current_active_user),
):
    """
    Create new gate
    """
    items = service.create(db, params=params)
    return items


@router.delete("/gates/{gate_id}", response_model=GateDto)
def delete_by_id(
    gate_id: int, current_user: UserModel = Depends(get_current_active_user), db: Session = Depends(get_db),
):
    """
    Delete gate by id
    """
    item = service.delete_by_id(db, id=gate_id)
    return item
