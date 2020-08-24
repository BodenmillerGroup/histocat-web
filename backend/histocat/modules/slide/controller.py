import os

from fastapi import APIRouter, Depends
from sqlalchemy.orm import Session
from starlette.responses import FileResponse

from histocat.api.db import get_db
from histocat.api.security import get_active_user
from histocat.modules.user.models import UserModel

from . import service
from .dto import SlideDto

router = APIRouter()


@router.get("/slides/{id}/image", responses={200: {"content": {"image/png": {}}}})
async def read_slide_image(
    id: int,
    # user: User = Depends(get_current_active_user),
    db: Session = Depends(get_db),
):
    """
    Get slide image by id
    """
    item = service.get(db, id=id)
    return FileResponse(
        os.path.join(item.location, "origin", f"{item.name}_s{item.origin_id}_slide.png"), media_type="image/png",
    )


@router.delete("/slides/{id}", response_model=SlideDto)
def delete_by_id(
    id: int, user: UserModel = Depends(get_active_user), db: Session = Depends(get_db),
):
    """
    Delete a specific slide by id
    """
    item = service.remove(db, id=id)
    return item