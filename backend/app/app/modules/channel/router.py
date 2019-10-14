import logging
from typing import List

import cv2
import numpy as np
import ujson
from fastapi import APIRouter, Depends
from imctools.io.ometiffparser import OmetiffParser
from sqlalchemy.orm import Session
from starlette.requests import Request
from starlette.responses import StreamingResponse, UJSONResponse

from app.api.utils.db import get_db
from app.api.utils.security import get_current_active_user
from app.core.image import (
    scale_image,
    colorize,
    apply_filter,
    draw_scalebar,
    draw_legend,
    draw_mask,
)
from app.core.redis_manager import redis_manager
from app.core.utils import stream_bytes
from app.modules.analysis.router import get_additive_image
from app.modules.user.db import User
from . import crud
from .models import ChannelModel, ChannelStatsModel, ChannelStackModel

logger = logging.getLogger(__name__)
router = APIRouter()


@router.get("/", response_model=List[ChannelModel])
def read_channels(
    db: Session = Depends(get_db),
    skip: int = 0,
    limit: int = 100,
    current_user: User = Depends(get_current_active_user),
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
    content = await redis_manager.cache.get(request.url.path)
    if content:
        return UJSONResponse(content=ujson.loads(content))

    item = crud.get(db, id=id)

    parser = OmetiffParser(item.acquisition.location)
    acq = parser.get_imc_acquisition()
    data = acq.get_img_by_metal(item.metal)

    hist, edges = np.histogram(data.ravel(), bins=bins)
    content = {"hist": hist.tolist(), "edges": edges.tolist()}
    await redis_manager.cache.set(request.url.path, ujson.dumps(content))
    return UJSONResponse(content=content)


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

    parser = OmetiffParser(item.acquisition.location)
    acq = parser.get_imc_acquisition()
    data = acq.get_img_by_metal(item.metal)

    levels = (
        (min, max)
        if min is not None and max is not None
        else (item.min_intensity, item.max_intensity)
    )
    data = scale_image(data, levels)

    color = f"#{color}" if color else "#ffffff"
    image = colorize(data, color)

    image = cv2.cvtColor(image.astype(data.dtype), cv2.COLOR_BGR2RGB)

    status, result = cv2.imencode(".png", image)
    return StreamingResponse(stream_bytes(result), media_type="image/png")


@router.post("/stack")
async def download_channel_stack(
    params: ChannelStackModel,
    current_user: User = Depends(get_current_active_user),
    db: Session = Depends(get_db),
):
    """
    Download channel stack (additive) image
    """
    additive_image, legend_labels = get_additive_image(db, params.channels)

    # TODO: Bright-field effect
    # additive_image = additive_image[..., ::-1]

    additive_image = cv2.cvtColor(additive_image, cv2.COLOR_BGR2RGB)

    if params.filter.apply:
        additive_image = apply_filter(additive_image, params.filter)

    if params.datasetId and params.mask and params.mask.apply:
        additive_image = draw_mask(additive_image, params.mask)

    if params.legend.apply:
        additive_image = draw_legend(additive_image, legend_labels, params.legend)

    if params.scalebar.apply:
        additive_image = draw_scalebar(additive_image, params.scalebar)

    format = params.format if params.format is not None else "png"
    status, result = cv2.imencode(
        f".{format}", additive_image.astype(int) if format == "tiff" else additive_image
    )
    return StreamingResponse(stream_bytes(result), media_type=f"image/{format}")
