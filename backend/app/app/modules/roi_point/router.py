from typing import List

from fastapi import APIRouter, Depends
from sqlalchemy.orm import Session

from app.api.utils.db import get_db
from app.api.utils.security import get_current_active_superuser, get_current_active_user
from app.modules.user.db import User
from . import crud
from .models import ROIPointModel

router = APIRouter()


@router.get("/", response_model=List[ROIPointModel])
def read_all(
    db: Session = Depends(get_db),
    skip: int = 0,
    limit: int = 100,
    current_user: User = Depends(get_current_active_superuser),
):
    """
    Retrieve ROIs
    """
    items = crud.get_multi(db, skip=skip, limit=limit)
    return items


@router.get("/{id}", response_model=ROIPointModel)
def read_by_id(
    id: int,
    current_user: User = Depends(get_current_active_user),
    db: Session = Depends(get_db),
):
    """
    Get a specific ROI by id
    """
    item = crud.get(db, id=id)
    return item
