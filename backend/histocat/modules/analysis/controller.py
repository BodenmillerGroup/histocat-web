import logging
import os
from typing import Optional, Sequence

import anndata as ad
import cv2
import numpy as np
import scanpy as sc
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
from histocat.modules.analysis.processors import phenograph, tsne, umap
from histocat.modules.dataset import service as dataset_service
from histocat.modules.gate import service as gate_service
from histocat.modules.result import service as result_service
from histocat.modules.user.models import UserModel

from ...io.dataset import ANNDATA_FILE_EXTENSION
from .dto import (
    PhenographDto,
    PhenographSubmissionDto,
    PlotSeriesDto,
    RegionChannelStatsDto,
    RegionStatsSubmissionDto,
    TsneDto,
    TsneSubmissionDto,
    UmapDto,
    UmapSubmissionDto,
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


@router.get("/analysis/scatterplot")
async def get_scatter_plot_data(
    dataset_id: int,
    marker_x: str,
    marker_y: str,
    result_id: Optional[int] = None,
    heatmap_type: Optional[str] = None,
    heatmap: Optional[str] = None,
    user: UserModel = Depends(get_active_user),
    db: Session = Depends(get_db),
):
    """Read scatter plot data from the dataset."""

    if result_id:
        result = result_service.get(db, id=result_id)
        if not result:
            raise HTTPException(status_code=404, detail="Result not found.")
        location = os.path.join(result.location, f"output{ANNDATA_FILE_EXTENSION}")
    else:
        dataset = dataset_service.get(db, id=dataset_id)
        if not dataset:
            raise HTTPException(status_code=HTTP_404_NOT_FOUND, detail="Dataset not found.")
        cell_input = dataset.meta.get("cell")
        location = cell_input.get("location")

    adata = ad.read_h5ad(location)

    output = {
        "acquisitionIds": adata.obs["AcquisitionId"].tolist(),
        "cellIds": adata.obs["CellId"].tolist(),
        "objectNumbers": adata.obs["ObjectNumber"].tolist(),
        "x": {"label": marker_x, "data": adata.X[:, adata.var.index == marker_x][:, 0].tolist(),},
        "y": {"label": marker_y, "data": adata.X[:, adata.var.index == marker_y][:, 0].tolist(),},
    }

    if heatmap:
        if heatmap_type == "channel":
            heatmap_values = adata.X[:, adata.var.index == heatmap]
            output["heatmap"] = {"label": heatmap, "heatmap_type": heatmap_type, "data": heatmap_values[:, 0].tolist()}
        elif heatmap_type == "neighbor" or heatmap_type == "clustering":
            heatmap_values = sc.get.obs_df(adata, keys=[heatmap])
            output["heatmap"] = {
                "label": heatmap,
                "heatmap_type": heatmap_type,
                "data": heatmap_values[heatmap].tolist(),
            }

    return ORJSONResponse(output)


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


@router.post("/analysis/tsne")
def submit_tsne(
    params: TsneSubmissionDto, user: UserModel = Depends(get_active_user), db: Session = Depends(get_db),
):
    """
    Start t-SNE data processing
    """
    worker.process_tsne.send(
        params.dataset_id,
        params.acquisition_ids,
        params.n_components,
        params.perplexity,
        params.learning_rate,
        params.iterations,
        params.theta,
        params.init,
        params.markers,
    )
    return ORJSONResponse({"status": "submitted"})


@router.get("/groups/{group_id}/results/{result_id}/tsne", response_model=TsneDto)
async def get_tsne_data(
    group_id: int,
    result_id: int,
    heatmap_type: Optional[str] = None,
    heatmap: Optional[str] = None,
    user: UserModel = Depends(get_active_user),
    db: Session = Depends(get_db),
):
    """
    Read t-SNE result data
    """

    content = tsne.get_tsne_result(db, result_id, heatmap_type, heatmap)
    return ORJSONResponse(content)


@router.post("/analysis/umap")
def submit_umap(
    params: UmapSubmissionDto, user: UserModel = Depends(get_active_user), db: Session = Depends(get_db),
):
    """
    Start UMAP data processing
    """
    worker.process_umap.send(
        params.dataset_id,
        params.acquisition_ids,
        params.n_components,
        params.n_neighbors,
        params.metric,
        params.min_dist,
        params.markers,
    )
    return ORJSONResponse({"status": "submitted"})


@router.get("/groups/{group_id}/results/{result_id}/umap", response_model=UmapDto)
async def get_umap_data(
    group_id: int,
    result_id: int,
    heatmap_type: Optional[str] = None,
    heatmap: Optional[str] = None,
    user: UserModel = Depends(get_active_user),
    db: Session = Depends(get_db),
):
    """
    Read UMAP result data
    """

    content = umap.get_umap_result(db, result_id, heatmap_type, heatmap)
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
