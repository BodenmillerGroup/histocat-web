from typing import Sequence

from fastapi import APIRouter, Depends, HTTPException
from sqlalchemy.orm import Session
from starlette.status import HTTP_404_NOT_FOUND

from histocat.api.db import get_db
from histocat.api.security import get_active_member, get_active_user
from histocat.modules.user.models import UserModel

from ..member.models import MemberModel
from . import service
from .dto import GateCreateDto, GateDto, GateUpdateDto

router = APIRouter()


@router.get("/gates/{gate_id}", response_model=GateDto)
def get_by_id(
    gate_id: int, user: UserModel = Depends(get_active_user), db: Session = Depends(get_db),
):
    """
    Get gate by id
    """
    item = service.get_by_id(db, id=gate_id)
    return item


@router.get("/datasets/{dataset_id}/gates", response_model=Sequence[GateDto])
def get_dataset_gates(
    dataset_id: int, user: UserModel = Depends(get_active_user), db: Session = Depends(get_db),
):
    """
    Get all dataset gates
    """
    item = service.get_dataset_gates(db, dataset_id=dataset_id)
    return item


@router.post("/gates", response_model=GateDto)
def create(
    *, db: Session = Depends(get_db), params: GateCreateDto, user: UserModel = Depends(get_active_user),
):
    """
    Create new gate
    """
    items = service.create(db, params=params)
    return items


@router.patch("/groups/{group_id}/gates/{gate_id}", response_model=GateDto)
def update(
    group_id: int,
    gate_id: int,
    params: GateUpdateDto,
    member: MemberModel = Depends(get_active_member),
    db: Session = Depends(get_db),
):
    """
    Update pipeline
    """
    item = service.get_by_id(db, id=gate_id)
    if not item:
        raise HTTPException(status_code=HTTP_404_NOT_FOUND, detail=f"Gate id:{gate_id} not found")
    item = service.update(db, item=item, params=params)
    return item


@router.delete("/gates/{gate_id}", response_model=int)
def delete_by_id(
    gate_id: int, user: UserModel = Depends(get_active_user), db: Session = Depends(get_db),
):
    """
    Delete gate by id
    """
    service.delete_by_id(db, id=gate_id)
    return gate_id
