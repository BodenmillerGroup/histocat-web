import logging
from typing import Optional, Sequence

import anndata as ad
import cv2
import numpy as np
from fastapi import APIRouter, Depends, HTTPException, Query
from fastapi.responses import ORJSONResponse
from imctools.io.ometiff.ometiffparser import OmeTiffParser
from skimage.measure import regionprops
from sqlalchemy.orm import Session
from starlette.requests import Request
from starlette.status import HTTP_404_NOT_FOUND

from histocat import worker
from histocat.api.db import get_db
from histocat.api.security import get_active_user
from histocat.core.image import colorize, scale_image
from histocat.modules.acquisition import service as acquisition_crud
from histocat.modules.acquisition.dto import ChannelStackDto
from histocat.modules.analysis.processors import phenograph
from histocat.modules.dataset import service as dataset_service
from histocat.modules.gate import service as gate_service
from histocat.modules.user.models import UserModel

from .dto import (
    PhenographDto,
    PhenographSubmissionDto,
    PlotSeriesDto,
    RegionChannelStatsDto,
    RegionStatsSubmissionDto,
)

logger = logging.getLogger(__name__)

router = APIRouter()

RESULT_TYPE_ORIGIN = "origin"
RESULT_TYPE_MASK = "mask"


def get_additive_image(db: Session, params: ChannelStackDto):
    additive_image: Optional[np.ndarray] = None

    acquisition = acquisition_crud.get_by_id(db, id=params.acquisitionId)
    if not acquisition:
        raise HTTPException(status_code=HTTP_404_NOT_FOUND, detail="Acquisition not found.")
    parser = OmeTiffParser(acquisition.location)
    acq = parser.get_acquisition_data()

    for channel in params.channels:
        data = acq.get_image_by_name(channel.name)
        item = acquisition.channels.get(channel.name)
        levels = (
            (channel.min, channel.max)
            if channel.min is not None and channel.max is not None
            else (item.get("min_intensity"), item.get("max_intensity"))
        )
        data = scale_image(data, levels)

        color = channel.color if channel.color else "#ffffff"
        image = colorize(data, color)

        if additive_image is None:
            additive_image = np.zeros(shape=(data.shape[0], data.shape[1], 4), dtype=np.float32)
        additive_image += image
    return np.clip(additive_image, 0, 1, out=additive_image)


@router.get("/analysis/boxplot", response_model=Sequence[PlotSeriesDto])
async def get_box_plot_data(
    dataset_id: int,
    gate_id: Optional[int] = None,
    acquisition_ids: Sequence[int] = Query(None),
    markers: Sequence[str] = Query(None),
    user: UserModel = Depends(get_active_user),
    db: Session = Depends(get_db),
):
    """
    Read box plot data from the dataset
    """

    dataset = dataset_service.get(db, id=dataset_id)
    cell_input = dataset.meta.get("cell")

    if not cell_input or (not gate_id and len(acquisition_ids) == 0):
        raise HTTPException(status_code=400, detail="The dataset does not have a proper input.")

    adata = ad.read_h5ad(cell_input.get("location"))

    content = []

    if gate_id:
        gate = gate_service.get_by_id(db, id=gate_id)
        adata = adata[adata.obs["CellId"].isin(gate.cell_ids)]
    else:
        adata = adata[adata.obs["AcquisitionId"].isin(acquisition_ids)]

    for marker in markers:
        content.append(
            {"label": marker, "data": adata.layers["exprs"][:, adata.var.index == marker][:, 0].tolist(),}
        )

    # await redis_manager.cache.set(request.url.path, ujson.dumps(content))
    return ORJSONResponse(content)


@router.post("/analysis/phenograph")
def submit_phenograph(
    params: PhenographSubmissionDto, user: UserModel = Depends(get_active_user), db: Session = Depends(get_db),
):
    """
    Start PhenoGraph data processing
    """

    worker.process_phenograph.send(
        params.dataset_id,
        params.acquisition_ids,
        params.markers,
        params.clustering_algo,
        params.nearest_neighbors,
        params.jaccard,
        params.primary_metric,
        params.min_cluster_size,
    )
    return ORJSONResponse({"status": "submitted"})


@router.get("/groups/{group_id}/results/{result_id}/phenograph", response_model=PhenographDto)
async def get_phenograph_data(
    group_id: int, result_id: int, user: UserModel = Depends(get_active_user), db: Session = Depends(get_db),
):
    """
    Read PhenoGraph result data
    """

    content = phenograph.get_phenograph_result(db, result_id)
    return ORJSONResponse(content)


@router.post("/analysis/region/stats", response_model=Sequence[RegionChannelStatsDto])
async def calculate_region_stats(
    params: RegionStatsSubmissionDto,
    request: Request,
    user: UserModel = Depends(get_active_user),
    db: Session = Depends(get_db),
):
    """
    Calculate region's statistics
    """

    acquisition = acquisition_crud.get_by_id(db, params.acquisition_id)
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
