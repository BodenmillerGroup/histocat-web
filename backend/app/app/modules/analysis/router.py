import logging
import os

import cv2
import numpy as np
import redis
from fastapi import APIRouter, Depends
from sqlalchemy.orm import Session
from starlette.responses import StreamingResponse

from app.api.utils.db import get_db
from app.api.utils.security import get_current_active_user
from app.core.image import scale_image, colorize, apply_filter, draw_scalebar
from app.core.utils import stream_bytes
from app.modules.user.db import User
from app.modules.channel import crud
from .models import AnalysisModel

logger = logging.getLogger(__name__)
r = redis.Redis(host="redis")

router = APIRouter()


@router.post("/segmentation")
async def download_analysis_image(
    params: AnalysisModel,
    current_user: User = Depends(get_current_active_user),
    db: Session = Depends(get_db),
):
    """
    Download analysis image
    """
    additive_image: np.ndarray = None
    for channel in params.channels:
        item = crud.get(db, id=channel.id)
        data = np.load(os.path.join(item.location, "origin.npy"))

        levels = (channel.min, channel.max) if channel.min is not None and channel.max is not None else (
            item.min_intensity, item.max_intensity)
        data = scale_image(data, levels)

        color = channel.color if channel.color else '#ffffff'
        image = colorize(data, color)

        if additive_image is None:
            additive_image = np.zeros(shape=(data.shape[0], data.shape[1], 4), dtype=data.dtype)
        additive_image += image

    # TODO: Bright-field effect
    # additive_image = additive_image[..., ::-1]

    additive_image = cv2.cvtColor(additive_image, cv2.COLOR_BGR2RGB)

    if params.filter.apply:
        additive_image = apply_filter(additive_image, params.filter)

    if params.scalebar.apply:
        additive_image = draw_scalebar(additive_image, params.scalebar)

    format = params.format if params.format is not None else 'png'
    status, result = cv2.imencode(f".{format}", additive_image.astype(int) if format == 'tiff' else additive_image)
    return StreamingResponse(stream_bytes(result), media_type=f"image/{format}")
