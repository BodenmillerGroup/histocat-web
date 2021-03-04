import logging
from typing import Sequence

import cv2
import numpy as np
from fastapi import APIRouter, Depends
from fastapi.responses import ORJSONResponse
from imctools.io.ometiff.ometiffparser import OmeTiffParser
from skimage.measure import regionprops
from sqlalchemy.orm import Session
from starlette.requests import Request

from histocat.api.db import get_db
from histocat.api.security import get_active_member
from histocat.core.acquisition import service as acquisition_service
from histocat.core.analysis.dto import RegionChannelStatsDto, RegionStatsSubmissionDto
from histocat.core.member.models import MemberModel

logger = logging.getLogger(__name__)

router = APIRouter()


@router.post("/groups/{group_id}/analysis/region/stats", response_model=Sequence[RegionChannelStatsDto])
async def calculate_region_stats(
    params: RegionStatsSubmissionDto,
    request: Request,
    member: MemberModel = Depends(get_active_member),
    db: Session = Depends(get_db),
):
    """
    Calculate region's statistics
    """

    acquisition = acquisition_service.get_by_id(db, params.acquisition_id)
    parser = OmeTiffParser(acquisition.location)
    acq = parser.get_acquisition_data()
    mask = None
    contour = np.array(params.region_polygon).astype(int)
    content = []
    for metal in acq.channel_names:
        channel_img = acq.get_image_by_name(metal)
        if mask is None:
            mask = np.zeros(channel_img.shape, np.uint8)
            mask = cv2.drawContours(mask, [contour], 0, 255, -1)
        props = regionprops(mask, intensity_image=channel_img, cache=True, coordinates=None)
        props = props[0]
        content.append(
            {
                "metal": metal,
                "min": float("{0:.3f}".format(props.min_intensity)),
                "max": float("{0:.3f}".format(props.max_intensity)),
                "mean": float("{0:.3f}".format(props.mean_intensity)),
            }
        )

    return ORJSONResponse(content)
