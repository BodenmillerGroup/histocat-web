import logging
from typing import Sequence

from fastapi import APIRouter, Depends, HTTPException
from sqlalchemy.orm import Session

from histocat.api.db import get_db
from histocat.api.security import get_active_member, get_active_user
from histocat.modules.member.models import MemberModel
from histocat.modules.user.models import UserModel

from . import service
from .dto import ResultDto, ResultUpdateDto

logger = logging.getLogger(__name__)
router = APIRouter()


@router.get("/groups/{group_id}/datasets/{dataset_id}/results", response_model=Sequence[ResultDto])
def get_dataset_results(
    group_id: int, dataset_id: int, db: Session = Depends(get_db), member: MemberModel = Depends(get_active_member),
):
    """
    Retrieve own results for specified dataset
    """
    items = service.get_dataset_results(db, dataset_id=dataset_id)
    return items


@router.get("/groups/{group_id}/results/{result_id}", response_model=ResultDto)
def get_by_id(
    group_id: int, result_id: int, user: UserModel = Depends(get_active_user), db: Session = Depends(get_db),
):
    """
    Get a specific result by id
    """
    item = service.get(db, id=result_id)
    return item


@router.patch("/groups/{group_id}/results/{result_id}", response_model=ResultDto)
def update(
    group_id: int,
    result_id: int,
    params: ResultUpdateDto,
    member: MemberModel = Depends(get_active_member),
    db: Session = Depends(get_db),
):
    """
    Update result
    """
    item = service.get(db, id=result_id)
    if not item:
        raise HTTPException(
            status_code=404, detail="Result not found",
        )
    item = service.update(db, item=item, params=params)
    return item


@router.delete("/groups/{group_id}/results/{result_id}", response_model=ResultDto)
def delete_by_id(
    group_id: int, result_id: int, user: UserModel = Depends(get_active_user), db: Session = Depends(get_db),
):
    """
    Delete a specific dataset by id
    """
    item = service.remove(db, id=result_id)
    return item
