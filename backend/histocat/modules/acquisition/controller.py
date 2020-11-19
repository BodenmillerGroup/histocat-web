import logging
import os
from io import BytesIO
from typing import Optional

import cv2
import imageio
import numpy as np
import orjson
import scanpy as sc
from fastapi import APIRouter, Depends, HTTPException
from fastapi.responses import ORJSONResponse
from imctools.io.ometiff.ometiffparser import OmeTiffParser
from skimage.util import img_as_ubyte
from sqlalchemy.orm import Session
from starlette.requests import Request
from starlette.responses import StreamingResponse
from starlette.status import HTTP_404_NOT_FOUND

from histocat.api.db import get_db
from histocat.api.security import get_active_member, get_active_user
from histocat.core.image import (
    apply_filter,
    colorize,
    draw_mask,
    draw_scalebar,
    scale_image,
)
from histocat.core.redis_manager import redis_manager
from histocat.core.utils import stream_bytes
from histocat.io.dataset import ANNDATA_FILE_EXTENSION
from histocat.modules.acquisition import service as acquisition_service
from histocat.modules.acquisition.dto import (
    ChannelStackDto,
    ChannelStatsDto,
    ChannelUpdateDto,
)
from histocat.modules.analysis.controller import get_additive_image
from histocat.modules.member.models import MemberModel
from histocat.modules.project.dto import ProjectFullDto
from histocat.modules.result import service as result_service
from histocat.modules.user.models import UserModel

logger = logging.getLogger(__name__)
router = APIRouter()


@router.get("/acquisitions/{acquisition_id}/{channel_name}/stats", response_model=ChannelStatsDto)
async def read_channel_stats(
    acquisition_id: int,
    channel_name: str,
    request: Request,
    bins: int = 40,
    user: UserModel = Depends(get_active_user),
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


@router.put("/groups/{group_id}/acquisitions/{acquisition_id}", response_model=ProjectFullDto)
def update(
    group_id: int,
    acquisition_id: int,
    params: ChannelUpdateDto,
    member: MemberModel = Depends(get_active_member),
    db: Session = Depends(get_db),
):
    """
    Update channel's custom label
    """
    item = acquisition_service.get_by_id(db, acquisition_id)
    if not item:
        raise HTTPException(
            status_code=HTTP_404_NOT_FOUND, detail="The acquisition with this id does not exist.",
        )
    item = acquisition_service.update_custom_label(db, item=item, params=params)
    return item


@router.get("/acquisitions/{acquisition_id}/{channel_name}/image", responses={200: {"content": {"image/png": {}}}})
async def read_channel_image(
    acquisition_id: int,
    channel_name: str,
    color: Optional[str] = None,
    min: Optional[float] = None,
    max: Optional[float] = None,
    # user: UserModel = Depends(get_current_active_user),
    db: Session = Depends(get_db),
):
    """Get channel image by name."""
    acquisition = acquisition_service.get_by_id(db, acquisition_id)

    parser = OmeTiffParser(acquisition.location)
    acq = parser.get_acquisition_data()
    data = acq.get_image_by_name(channel_name)

    channel = acquisition.channels[channel_name]

    levels = (
        (min, max)
        if min is not None and max is not None
        else (channel.get("min_intensity"), channel.get("max_intensity"))
    )
    data = scale_image(data, levels)

    color = f"#{color}" if color else "#ffffff"
    image = colorize(data, color)

    image = img_as_ubyte(image)
    image = cv2.cvtColor(image, cv2.COLOR_BGR2RGB)
    status, buffer = cv2.imencode(".png", image)
    return StreamingResponse(stream_bytes(buffer), media_type="image/png")


@router.post("/acquisitions/stack")
async def download_channel_stack(
    params: ChannelStackDto, user: UserModel = Depends(get_active_user), db: Session = Depends(get_db),
):
    """Download channel stack (additive) image."""
    additive_image = get_additive_image(db, params)

    # TODO: Bright-field effect
    # additive_image = additive_image[..., ::-1]

    if params.filter.apply:
        additive_image = apply_filter(additive_image, params.filter)

    if params.datasetId and params.mask and params.mask.apply:
        heatmap_dict = None
        if params.mask.resultId and params.mask.colorsType and params.mask.colorsName:
            result = result_service.get(db, id=params.mask.resultId)
            location = os.path.join(result.location, f"output{ANNDATA_FILE_EXTENSION}")
            adata = sc.read_h5ad(location)
            adata = adata[adata.obs["AcquisitionId"] == params.acquisitionId]

            heatmap_values = None
            if params.mask.colorsType == "marker":
                heatmap_values = adata.X[:, adata.var.index == params.mask.colorsName]
                heatmap_dict = dict(zip(adata.obs["ObjectNumber"], heatmap_values[:, 0].tolist()))
                heatmap_dict.pop('0', None)
            elif params.mask.colorsType == "neighbor" or params.mask.colorsType == "clustering":
                heatmap_values = sc.get.obs_df(adata, keys=[params.mask.colorsName])
                heatmap_dict = dict(zip(adata.obs["ObjectNumber"], heatmap_values[params.mask.colorsName].astype("uint")))
                heatmap_dict.pop('0', None)

        additive_image = draw_mask(additive_image, params.mask, heatmap_dict)

    if params.scalebar.apply:
        additive_image = draw_scalebar(additive_image, params.scalebar)

    format = params.format if params.format is not None else "png"

    # buffer = BytesIO()
    # imageio.imwrite(buffer, additive_image, format=format)
    # buffer.seek(0)
    #
    # return StreamingResponse(buffer, media_type=f"image/{format}")

    additive_image = img_as_ubyte(additive_image)
    additive_image = cv2.cvtColor(additive_image, cv2.COLOR_BGR2RGB)
    status, buffer = cv2.imencode(f".{format}", additive_image.astype(int) if format == "tiff" else additive_image)
    return StreamingResponse(stream_bytes(buffer), media_type=f"image/{format}")
