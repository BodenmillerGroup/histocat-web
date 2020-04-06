import os
from typing import List

from fastapi import APIRouter, Depends
from sqlalchemy.orm import Session
from starlette.responses import FileResponse

from app.api.utils.db import get_db
from app.api.utils.security import get_current_active_user
from app.modules.user.models import User

from . import service
from .dto import SlideDto

router = APIRouter()


@router.get("/", response_model=List[SlideDto])
def read_all(
    db: Session = Depends(get_db),
    skip: int = 0,
    limit: int = 100,
    current_user: User = Depends(get_current_active_user),
):
    """
    Retrieve slides
    """
    items = service.get_multi(db, skip=skip, limit=limit)
    return items


@router.get("/{id}", response_model=SlideDto)
def read_by_id(
    id: int, current_user: User = Depends(get_current_active_user), db: Session = Depends(get_db),
):
    """
    Get a specific slide by id
    """
    item = service.get(db, id=id)
    return item


@router.get("/{id}/image", responses={200: {"content": {"image/png": {}}}})
async def read_slide_image(
    id: int,
    # current_user: User = Depends(get_current_active_user),
    db: Session = Depends(get_db),
):
    """
    Get slide image by id
    """
    item = service.get(db, id=id)
    return FileResponse(
        os.path.join(item.location, "origin", f"{item.name}_s{item.origin_id}_slide.png"), media_type="image/png",
    )


@router.delete("/{id}", response_model=SlideDto)
def delete_by_id(
    id: int, current_user: User = Depends(get_current_active_user), db: Session = Depends(get_db),
):
    """
    Delete a specific slide by id
    """
    item = service.remove(db, id=id)
    return item
