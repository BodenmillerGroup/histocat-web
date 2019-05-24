from typing import List

from fastapi import APIRouter, Depends, HTTPException
from sqlalchemy.orm import Session

from app.api.utils.db import get_db
from app.api.utils.security import get_current_active_superuser, get_current_active_user
from app.modules.user.crud import is_superuser
from app.modules.user.db import User

from . import crud
from .models import AcquisitionCreateModel, AcquisitionModel, AcquisitionUpdateModel

router = APIRouter()


@router.get("/", response_model=List[AcquisitionModel])
def read_acquisitions(
    db: Session = Depends(get_db),
    skip: int = 0,
    limit: int = 100,
    current_user: User = Depends(get_current_active_superuser),
):
    """
    Retrieve acquisitions
    """
    items = crud.get_multi(db, skip=skip, limit=limit)
    return items


@router.post("/", response_model=AcquisitionModel)
def create_acquisition(
    *,
    db: Session = Depends(get_db),
    params: AcquisitionCreateModel,
    current_user: User = Depends(get_current_active_superuser),
):
    """
    Create new acquisition
    """
    item = crud.get_by_name(db, name=params.name)
    if item:
        raise HTTPException(
            status_code=400,
            detail="The acquisition with this name already exists in the system.",
        )
    item = crud.create(db, params=params)
    return item


@router.get("/{id}", response_model=AcquisitionModel)
def read_acquisition_by_id(
    id: int,
    current_user: User = Depends(get_current_active_user),
    db: Session = Depends(get_db),
):
    """
    Get a specific acquisition by id
    """
    item = crud.get(db, id=id)
    return item


@router.put("/{id}", response_model=AcquisitionModel)
def update_acquisition(
    *,
    db: Session = Depends(get_db),
    id: int,
    params: AcquisitionUpdateModel,
    current_user: User = Depends(get_current_active_superuser),
):
    """
    Update a slide
    """
    item = crud.get(db, id=id)

    if not item:
        raise HTTPException(status_code=404, detail="Item not found")
    if not is_superuser(current_user):
        raise HTTPException(status_code=400, detail="Not enough permissions")
    item = crud.update(db, item=item, params=params)
    return item
