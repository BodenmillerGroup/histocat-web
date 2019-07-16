import os
from io import BytesIO
from typing import List

import cv2
import numpy as np
import redis
import ujson
from fastapi import APIRouter, Depends
from sqlalchemy.orm import Session
from starlette.requests import Request
from starlette.responses import StreamingResponse, UJSONResponse

from app.api.utils.db import get_db
from app.api.utils.security import get_current_active_superuser, get_current_active_user
from app.core.utils import Color, colorize, scale_image
from app.modules.user.db import User

from . import crud
from .models import ChannelModel, ChannelStatsModel

r = redis.Redis(host="redis")

router = APIRouter()


@router.get("/", response_model=List[ChannelModel])
def read_channels(
    db: Session = Depends(get_db),
    skip: int = 0,
    limit: int = 100,
    current_user: User = Depends(get_current_active_superuser),
):
    """
    Retrieve channels
    """
    items = crud.get_multi(db, skip=skip, limit=limit)
    return items


@router.get("/{id}", response_model=ChannelModel)
def read_channel_by_id(
    id: int,
    current_user: User = Depends(get_current_active_user),
    db: Session = Depends(get_db),
):
    """
    Get a specific channel by id
    """
    item = crud.get(db, id=id)
    return item


async def stream_image(record: bytes, chunk_size: int = 65536):
    with BytesIO(record) as stream:
        data = stream.read(chunk_size)
        while data:
            yield data
            data = stream.read(chunk_size)


@router.get("/{id}/image", responses={200: {"content": {"image/png": {}}}})
async def read_channel_image(
    id: int,
    color: str = None,
    min: float = None,
    max: float = None,
    # current_user: User = Depends(get_current_active_user),
    db: Session = Depends(get_db),
):
    """
    Get channel image by id
    """
    item = crud.get(db, id=id)

    data = np.load(os.path.join(item.location, "origin.npy"))
    if min is not None and max is not None:
        data = scale_image(data, (item.min_intensity, item.max_intensity), (min, max))
    clr = Color[color] if color else None
    img = colorize(data, clr) if clr else data
    png = cv2.imencode(".png", img)[1]
    return StreamingResponse(
        stream_image(png), media_type="image/png", headers={"Cache-Control": "private"}
    )


@router.get("/{id}/stats", response_model=ChannelStatsModel)
async def read_channel_stats(
    id: int,
    request: Request,
    bins: int = 100,
    current_user: User = Depends(get_current_active_user),
    db: Session = Depends(get_db),
):
    """
    Get channel stats by id
    """
    content = r.get(request.url.path)
    if content:
        return UJSONResponse(
            content=ujson.loads(content), headers={"Cache-Control": "private"}
        )

    item = crud.get(db, id=id)

    data = np.load(os.path.join(item.location, "origin.npy"))
    hist, edges = np.histogram(data.ravel(), bins=bins)
    content = {"hist": hist.tolist(), "edges": edges.tolist()}
    r.set(request.url.path, ujson.dumps(content))
    return UJSONResponse(content=content, headers={"Cache-Control": "private"})
