import logging
from typing import List, Optional, Tuple

import cv2
import numpy as np
from fastapi import APIRouter, Depends
from imctools.io.ometiffparser import OmetiffParser
from matplotlib.colors import to_rgba
from sqlalchemy.orm import Session
from starlette.requests import Request
from starlette.responses import StreamingResponse, UJSONResponse

from app.api.utils.db import get_db
from app.api.utils.security import get_current_active_user
from app.core.image import scale_image, colorize, apply_filter, draw_scalebar, get_mask, apply_morphology
from app.core.utils import stream_bytes
from app.modules.channel import crud
from app.modules.channel.models import ChannelSettingsModel
from app.modules.user.db import User
from .models import AnalysisModel

logger = logging.getLogger(__name__)

router = APIRouter()

RESULT_TYPE_ORIGIN = 'origin'
RESULT_TYPE_MASK = 'mask'


def get_additive_image(db: Session, channels: List[ChannelSettingsModel]):
    additive_image: Optional[np.ndarray] = None
    legend_labels: List[Tuple[str, str, float]] = list()

    item = channels[0]
    first = crud.get(db, id=item.id)
    parser = OmetiffParser(first.location)
    acq = parser.get_imc_acquisition()

    for channel in channels:
        item = crud.get(db, id=channel.id)

        data = acq.get_img_by_metal(item.metal)

        levels = (channel.min, channel.max) if channel.min is not None and channel.max is not None else (
            item.min_intensity, item.max_intensity)
        data = scale_image(data, levels)

        color = channel.color if channel.color else '#ffffff'
        image = colorize(data, color)

        label = channel.customLabel if channel.customLabel else item.label
        legend_labels.append((label, color, levels[1]))

        if additive_image is None:
            additive_image = np.zeros(shape=(data.shape[0], data.shape[1], 4), dtype=data.dtype)
        additive_image += image
    return additive_image, legend_labels


@router.post("/segmentation/image")
async def produce_segmentation_image(
    params: AnalysisModel,
    current_user: User = Depends(get_current_active_user),
    db: Session = Depends(get_db),
):
    """
    Produce segmentation image
    """
    additive_image, _ = get_additive_image(db, params.channels)

    if params.filter.apply:
        additive_image = apply_filter(additive_image, params.filter)

    if params.scalebar.apply:
        additive_image = draw_scalebar(additive_image, params.scalebar)

    mask = get_mask(additive_image, params.settings)

    if params.settings.iterations > 0:
        mask = apply_morphology(mask, params.settings)

    if params.settings.result_type == RESULT_TYPE_ORIGIN:
        r, g, b, a = to_rgba(params.settings.mask_color)
        additive_image[mask == 0] = (r * 255, g * 255, b * 255, a * 255)
        additive_image = cv2.cvtColor(additive_image, cv2.COLOR_BGRA2RGBA)
    elif params.settings.result_type == RESULT_TYPE_MASK:
        additive_image = mask
    else:
        pass

    format = params.format if params.format is not None else 'png'
    status, result = cv2.imencode(f".{format}", additive_image.astype(int) if format == 'tiff' else additive_image)
    return StreamingResponse(stream_bytes(result), media_type=f"image/{format}")


@router.post("/segmentation/contours")
async def produce_segmentation_contours(
    params: AnalysisModel,
    request: Request,
    current_user: User = Depends(get_current_active_user),
    db: Session = Depends(get_db),
):
    """
    Produce segmentation image
    """
    additive_image, _ = get_additive_image(db, params.channels)

    if params.filter.apply:
        additive_image = apply_filter(additive_image, params.filter)

    mask = get_mask(additive_image, params.settings)

    if params.settings.iterations > 0:
        mask = apply_morphology(mask, params.settings)

    contours0, hierarchy = cv2.findContours(cv2.flip(mask, 0), mode=cv2.RETR_LIST, method=cv2.CHAIN_APPROX_SIMPLE)
    contours = [cv2.approxPolyDP(cnt, 3, True).tolist() for cnt in contours0]
    return UJSONResponse(content=contours)
