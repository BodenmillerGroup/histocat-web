import os
from typing import List

from fastapi import APIRouter, Depends
from sqlalchemy.orm import Session
from starlette.responses import FileResponse

from app.api.utils.db import get_db
from app.api.utils.security import get_current_active_superuser, get_current_active_user
from app.modules.user.db import User
from . import crud
from .models import PanoramaModel

router = APIRouter()


@router.get("/", response_model=List[PanoramaModel])
def read_all(
    db: Session = Depends(get_db),
    skip: int = 0,
    limit: int = 100,
    current_user: User = Depends(get_current_active_superuser),
):
    """
    Retrieve panoramas
    """
    items = crud.get_multi(db, skip=skip, limit=limit)
    return items


@router.get("/{id}", response_model=PanoramaModel)
def read_by_id(
    id: int,
    current_user: User = Depends(get_current_active_user),
    db: Session = Depends(get_db),
):
    """
    Get a specific panorama by id
    """
    item = crud.get(db, id=id)
    return item


@router.get("/{id}/image", responses={200: {"content": {"image/png": {}}}})
async def read_panorama_image(
    id: int,
    # current_user: User = Depends(get_current_active_user),
    db: Session = Depends(get_db),
):
    """
    Get panorama image by id
    """
    item = crud.get(db, id=id)
    slide = item.slide
    return FileResponse(os.path.join(item.slide.location, "origin", f"{slide.name}_s{slide.origin_id}_p{item.origin_id}_pano.png"), media_type="image/png")
