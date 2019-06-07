from typing import List

from fastapi import APIRouter, Depends, HTTPException
from sqlalchemy.orm import Session

from app.api.utils.db import get_db
from app.api.utils.security import get_current_active_superuser, get_current_active_user
from app.modules.user.db import User

from . import crud
from .models import SlideCreateModel, SlideModel

router = APIRouter()


@router.get("/", response_model=List[SlideModel])
def read_slides(
    db: Session = Depends(get_db),
    skip: int = 0,
    limit: int = 100,
    current_user: User = Depends(get_current_active_superuser),
):
    """
    Retrieve slides
    """
    items = crud.get_multi(db, skip=skip, limit=limit)
    return items


@router.post("/", response_model=SlideModel)
def create_slide(
    *,
    db: Session = Depends(get_db),
    params: SlideCreateModel,
    current_user: User = Depends(get_current_active_superuser),
):
    """
    Create new slide
    """
    item = crud.get_by_name(db, name=params.name)
    if item:
        raise HTTPException(
            status_code=400,
            detail="The slide with this name already exists in the system.",
        )
    item = crud.create(db, params=params)
    return item


@router.get("/{id}", response_model=SlideModel)
def read_slide_by_id(
    id: int,
    current_user: User = Depends(get_current_active_user),
    db: Session = Depends(get_db),
):
    """
    Get a specific slide by id
    """
    item = crud.get(db, id=id)
    return item
