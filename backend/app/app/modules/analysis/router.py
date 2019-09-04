import logging
from typing import List, Optional, Tuple

import cv2
import numpy as np
import pandas as pd
import ujson
from fastapi import APIRouter, Depends, HTTPException, Query
from imctools.io.ometiffparser import OmetiffParser
from matplotlib.colors import to_rgba
from sqlalchemy.orm import Session
from starlette.requests import Request
from starlette.responses import StreamingResponse, UJSONResponse

from app.api.utils.db import get_db
from app.api.utils.security import get_current_active_user
from app.core.image import scale_image, colorize, apply_filter, draw_scalebar, get_mask, apply_morphology
from app.core.utils import stream_bytes
from app.modules.channel import crud as channel_crud
from app.modules.dataset import crud as dataset_crud
from app.modules.channel.models import ChannelSettingsModel
from app.modules.user.db import User
from .models import AnalysisModel, ScatterPlotModel, BoxPlotModel, PlotSeriesModel
from app.core.redis_manager import redis_manager

logger = logging.getLogger(__name__)

router = APIRouter()

RESULT_TYPE_ORIGIN = 'origin'
RESULT_TYPE_MASK = 'mask'


def get_additive_image(db: Session, channels: List[ChannelSettingsModel]):
    additive_image: Optional[np.ndarray] = None
    legend_labels: List[Tuple[str, str, float]] = list()

    item = channels[0]
    first = channel_crud.get(db, id=item.id)
    parser = OmetiffParser(first.acquisition.location)
    acq = parser.get_imc_acquisition()

    for channel in channels:
        item = channel_crud.get(db, id=channel.id)

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


@router.get("/scatterplot", response_model=ScatterPlotModel)
async def read_scatter_plot_data(
    request: Request,
    dataset_id: int,
    acquisition_id: int,
    marker_x: str,
    marker_y: str,
    marker_z: Optional[str] = None,
    marker_color: str = None,
    current_user: User = Depends(get_current_active_user),
    db: Session = Depends(get_db),
):
    """
    Read scatter plot data from the dataset
    """
    # content = await redis_manager.cache.get(request.url.path)
    # if content:
    #     return UJSONResponse(content=ujson.loads(content))

    dataset = dataset_crud.get(db, id=dataset_id)
    cell_artifact = dataset.artifacts.get("cell")
    channel_map = dataset.artifacts.get("channel_map")
    image_map = dataset.artifacts.get("image_map")
    image_number = image_map.get(str(acquisition_id))
    if not cell_artifact or not image_number or not channel_map:
        raise HTTPException(
            status_code=400,
            detail="The dataset does not have proper artifacts.",
        )

    content = {}
    df = pd.read_feather(cell_artifact.get("location"))
    df = df[df["ImageNumber"] == image_number]
    content["x"] = {
        "marker": marker_x,
        "data": df[f'Intensity_MeanIntensity_FullStack_c{channel_map[marker_x]}'] * 2**16
    }
    content["y"] = {
        "marker": marker_y,
        "data": df[f'Intensity_MeanIntensity_FullStack_c{channel_map[marker_y]}'] * 2**16
    }

    if marker_z:
        content["z"] = {
            "marker": marker_z,
            "data": df[f'Intensity_MeanIntensity_FullStack_c{channel_map[marker_z]}'] * 2**16
        }

    # await redis_manager.cache.set(request.url.path, ujson.dumps(content))
    return UJSONResponse(content=content)


@router.get("/boxplot", response_model=List[PlotSeriesModel])
async def read_box_plot_data(
    request: Request,
    dataset_id: int,
    acquisition_id: int,
    markers: List[str] = Query(None),
    current_user: User = Depends(get_current_active_user),
    db: Session = Depends(get_db),
):
    """
    Read box plot data from the dataset
    """
    # content = await redis_manager.cache.get(request.url.path)
    # if content:
    #     return UJSONResponse(content=ujson.loads(content))

    dataset = dataset_crud.get(db, id=dataset_id)
    cell_artifact = dataset.artifacts.get("cell")
    channel_map = dataset.artifacts.get("channel_map")
    image_map = dataset.artifacts.get("image_map")
    image_number = image_map.get(str(acquisition_id))
    if not cell_artifact or not image_number or not channel_map:
        raise HTTPException(
            status_code=400,
            detail="The dataset does not have proper artifacts.",
        )

    content = []
    df = pd.read_feather(cell_artifact.get("location"))
    df = df[df["ImageNumber"] == image_number]

    for marker in markers:
        content.append({
            "marker": marker,
            "data": df[f'Intensity_MeanIntensity_FullStack_c{channel_map[marker]}'] * 2 ** 16
        })

    # await redis_manager.cache.set(request.url.path, ujson.dumps(content))
    return UJSONResponse(content=content)
