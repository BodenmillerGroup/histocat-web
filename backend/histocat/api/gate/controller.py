from typing import Sequence

from fastapi import APIRouter, Depends, HTTPException
from sqlalchemy.orm import Session
from starlette.status import HTTP_404_NOT_FOUND

from histocat.api.db import get_db
from histocat.api.security import get_active_member
from histocat.core.gate import service
from histocat.core.gate.dto import GateCreateDto, GateDto, GateUpdateDto
from histocat.core.member.models import MemberModel

router = APIRouter()


@router.get("/groups/{group_id}/gates/{gate_id}", response_model=GateDto)
def get_by_id(
    group_id: int, gate_id: int, member: MemberModel = Depends(get_active_member), db: Session = Depends(get_db),
):
    """
    Get gate by id
    """
    item = service.get_by_id(db, id=gate_id)
    return item


@router.get("/groups/{group_id}/datasets/{dataset_id}/gates", response_model=Sequence[GateDto])
def get_dataset_gates(
    group_id: int, dataset_id: int, member: MemberModel = Depends(get_active_member), db: Session = Depends(get_db),
):
    """
    Get all dataset gates
    """
    item = service.get_dataset_gates(db, dataset_id=dataset_id)
    return item


@router.post("/groups/{group_id}/gates", response_model=GateDto)
def create(
    group_id: int,
    params: GateCreateDto,
    db: Session = Depends(get_db),
    member: MemberModel = Depends(get_active_member),
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
    Update gate
    """
    item = service.get_by_id(db, id=gate_id)
    if not item:
        raise HTTPException(status_code=HTTP_404_NOT_FOUND, detail=f"Gate id:{gate_id} not found")
    item = service.update(db, item=item, params=params)
    return item


@router.delete("/groups/{group_id}/gates/{gate_id}", response_model=int)
def delete_by_id(
    group_id: int, gate_id: int, member: MemberModel = Depends(get_active_member), db: Session = Depends(get_db),
):
    """
    Delete gate by id
    """
    service.delete_by_id(db, id=gate_id)
    return gate_id
