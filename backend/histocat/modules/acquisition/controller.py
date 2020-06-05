import logging
from typing import Optional

import cv2
import numpy as np
import orjson
from fastapi import APIRouter, Depends
from fastapi.responses import ORJSONResponse
from imctools.io.ometiff.ometiffparser import OmeTiffParser
from sqlalchemy.orm import Session
from starlette.requests import Request
from starlette.responses import StreamingResponse

from histocat.api.utils.db import get_db
from histocat.api.utils.security import get_current_active_user
from histocat.core.image import (
    apply_filter,
    colorize,
    draw_mask,
    draw_scalebar,
    scale_image,
)
from histocat.core.redis_manager import redis_manager
from histocat.core.utils import stream_bytes
from histocat.modules.acquisition import service as acquisition_service
from histocat.modules.acquisition.dto import ChannelStackDto, ChannelStatsDto
from histocat.modules.analysis.controller import get_additive_image
from histocat.modules.user.models import UserModel

logger = logging.getLogger(__name__)
router = APIRouter()


@router.get("/{acquisition_id}/{channel_name}/stats", response_model=ChannelStatsDto)
async def read_channel_stats(
    acquisition_id: int,
    channel_name: str,
    request: Request,
    bins: int = 40,
    current_user: UserModel = Depends(get_current_active_user),
    db: Session = Depends(get_db),
):
    """Get channel stats by name."""
    content = await redis_manager.cache.get(request.url.path)
    if content:
        return orjson.loads(content)

    acquisition = acquisition_service.get_by_id(db, acquisition_id)

    parser = OmeTiffParser(acquisition.location)
    acq = parser.get_acquisition_data()
    data = acq.get_image_by_name(channel_name)

    # TODO: check if the transformation is really needed
    # data = np.arcsinh(data / 5, out=data)

    hist, _ = np.histogram(data.ravel(), bins=bins)
    content = {"bins": hist.tolist()}
    await redis_manager.cache.set(request.url.path, orjson.dumps(content))
    return ORJSONResponse(content)


@router.get("/{acquisition_id}/{channel_name}/image", responses={200: {"content": {"image/png": {}}}})
async def read_channel_image(
    acquisition_id: int,
    channel_name: str,
    color: Optional[str] = None,
    min: Optional[float] = None,
    max: Optional[float] = None,
    # current_user: UserModel = Depends(get_current_active_user),
    db: Session = Depends(get_db),
):
    """Get channel image by name."""
    acquisition = acquisition_service.get_by_id(db, acquisition_id)

    parser = OmeTiffParser(acquisition.location)
    acq = parser.get_acquisition_data()
    data = acq.get_image_by_name(channel_name)

    channel = acquisition.channels[channel_name]

    # lmin = float(data.min())
    # lmax = float(data.max())
    # data = np.floor((data - lmin) / (lmax - lmin) * 255.)

    levels = (
        (min, max)
        if min is not None and max is not None
        else (channel.get("min_intensity"), channel.get("max_intensity"))
    )
    data = scale_image(data, levels)

    color = f"#{color}" if color else "#ffffff"
    image = colorize(data, color)

    image = cv2.cvtColor(image.astype(data.dtype), cv2.COLOR_BGR2RGB)

    status, result = cv2.imencode(".png", image)
    return StreamingResponse(stream_bytes(result), media_type="image/png")


@router.post("/stack")
async def download_channel_stack(
    params: ChannelStackDto, current_user: UserModel = Depends(get_current_active_user), db: Session = Depends(get_db),
):
    """Download channel stack (additive) image."""
    additive_image, legend_labels = get_additive_image(db, params)

    # TODO: Bright-field effect
    # additive_image = additive_image[..., ::-1]

    additive_image = cv2.cvtColor(additive_image, cv2.COLOR_BGR2RGB)

    if params.filter.apply:
        additive_image = apply_filter(additive_image, params.filter)

    if params.datasetId and params.mask and params.mask.apply:
        additive_image = draw_mask(additive_image, params.mask)

    if params.scalebar.apply:
        additive_image = draw_scalebar(additive_image, params.scalebar)

    format = params.format if params.format is not None else "png"
    status, result = cv2.imencode(f".{format}", additive_image.astype(int) if format == "tiff" else additive_image)
    return StreamingResponse(stream_bytes(result), media_type=f"image/{format}")
