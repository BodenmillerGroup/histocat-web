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
from app.core.utils import colorize, scale_image, apply_filter
from app.modules.user.db import User

from . import crud
from .models import ChannelModel, ChannelStatsModel, ChannelStackModel

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

    levels = (min, max) if min is not None and max is not None else (item.min_intensity, item.max_intensity)
    data = scale_image(data, levels)

    color = f'#{color}' if color else None
    image = colorize(data, color)

    image = cv2.cvtColor(image.astype(data.dtype), cv2.COLOR_BGR2RGB)

    status, result = cv2.imencode(".png", image)
    return StreamingResponse(stream_image(result), media_type="image/png")


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
        return UJSONResponse(content=ujson.loads(content))

    item = crud.get(db, id=id)

    data = np.load(os.path.join(item.location, "origin.npy"))
    hist, edges = np.histogram(data.ravel(), bins=bins)
    content = {"hist": hist.tolist(), "edges": edges.tolist()}
    r.set(request.url.path, ujson.dumps(content))
    return UJSONResponse(content=content)


@router.post("/stack")
async def download_channel_stack(
    params: ChannelStackModel,
    current_user: User = Depends(get_current_active_user),
    db: Session = Depends(get_db),
):
    """
    Download channel stack (additive) image
    """
    additive_image: np.ndarray = None
    for channel in params.channels:
        item = crud.get(db, id=channel.id)
        data = np.load(os.path.join(item.location, "origin.npy"))

        levels = (channel.min, channel.max) if channel.min is not None and channel.max is not None else (item.min_intensity, item.max_intensity)
        data = scale_image(data, levels)

        color = channel.color if channel.color else None
        image = colorize(data, color)

        if additive_image is None:
            additive_image = np.zeros(shape=(data.shape[0], data.shape[1], 4), dtype=data.dtype)
        additive_image += image

    # TODO: Bright-field effect
    # additive_image = additive_image[..., ::-1]

    additive_image = cv2.cvtColor(additive_image, cv2.COLOR_BGR2RGB)

    if params.filter.apply:
        additive_image = apply_filter(additive_image, params.filter)

    format = params.format if params.format is not None else 'png'
    status, result = cv2.imencode(f".{format}", additive_image.astype(int) if format == 'tiff' else additive_image)
    return StreamingResponse(stream_image(result), media_type=f"image/{format}")
