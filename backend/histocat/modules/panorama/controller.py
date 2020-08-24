import os

from fastapi import APIRouter, Depends
from sqlalchemy.orm import Session
from starlette.responses import FileResponse

from histocat.api.db import get_db

from . import service

router = APIRouter()


@router.get("/panoramas/{id}/image", responses={200: {"content": {"image/png": {}}}})
async def read_panorama_image(
    id: int,
    # user: User = Depends(get_current_active_user),
    db: Session = Depends(get_db),
):
    """
    Get panorama image by id
    """
    item = service.get(db, id=id)
    slide = item.slide
    return FileResponse(
        os.path.join(item.slide.location, "origin", f"{slide.name}_s{slide.origin_id}_p{item.origin_id}_pano.png",),
        media_type="image/png",
    )
