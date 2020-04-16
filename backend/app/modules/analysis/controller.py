import logging
from typing import List, Optional, Tuple

import cv2
import numpy as np
import pandas as pd
from fastapi import APIRouter, Depends, HTTPException, Query
from fastapi.responses import ORJSONResponse
from imctools.io.ometiff.ometiffparser import OmeTiffParser
from matplotlib.colors import to_rgba
from skimage.measure import regionprops
from sqlalchemy.orm import Session
from starlette.requests import Request
from starlette.responses import StreamingResponse

from app import worker
from app.api.utils.db import get_db
from app.api.utils.security import get_current_active_user
from app.core.image import (
    apply_filter,
    apply_morphology,
    colorize,
    draw_scalebar,
    get_mask,
    scale_image,
)
from app.core.utils import stream_bytes
from app.modules.acquisition import service as acquisition_crud
from app.modules.acquisition.dto import ChannelStackDto
from app.modules.analysis.processors import pca, phenograph, tsne, umap
from app.modules.dataset import service as dataset_crud
from app.modules.user.models import UserModel

from .dto import (
    AnalysisDto,
    PcaDto,
    PhenographDto,
    PhenographSubmissionDto,
    PlotSeriesDto,
    RegionChannelStatsDto,
    RegionStatsSubmissionDto,
    ScatterPlotDto,
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
    legend_labels: List[Tuple[str, str, float]] = list()

    acquisition = acquisition_crud.get_by_id(db, id=params.acquisitionId)
    if not acquisition:
        raise HTTPException(status_code=400, detail="Acquisition not found.")
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

        label = channel.customLabel if channel.customLabel else item.label
        legend_labels.append((label, color, levels[1]))

        if additive_image is None:
            additive_image = np.zeros(shape=(data.shape[0], data.shape[1], 4), dtype=data.dtype)
        additive_image += image
    return additive_image, legend_labels


@router.post("/segmentation/image")
async def produce_segmentation_image(
    params: AnalysisDto, current_user: UserModel = Depends(get_current_active_user), db: Session = Depends(get_db),
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

    format = params.format if params.format is not None else "png"
    status, result = cv2.imencode(f".{format}", additive_image.astype(int) if format == "tiff" else additive_image)
    return StreamingResponse(stream_bytes(result), media_type=f"image/{format}")


@router.post("/segmentation/contours")
async def produce_segmentation_contours(
    params: AnalysisDto,
    request: Request,
    current_user: UserModel = Depends(get_current_active_user),
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
    return ORJSONResponse(contours)


@router.get("/scatterplot")
async def read_scatter_plot_data(
    dataset_id: int,
    marker_x: str,
    marker_y: str,
    acquisition_ids: List[int] = Query(None),
    marker_z: Optional[str] = None,
    heatmap_type: Optional[str] = None,
    heatmap: Optional[str] = None,
    current_user: UserModel = Depends(get_current_active_user),
    db: Session = Depends(get_db),
):
    """
    Read scatter plot data from the dataset
    """
    # content = await redis_manager.cache.get(request.url.path)
    # if content:
    #     return UJSONResponse(content=ujson.loads(content))

    dataset = dataset_crud.get(db, id=dataset_id)
    cell_input = dataset.input.get("cell")
    channel_map = dataset.input.get("channel_map")
    image_map = dataset.input.get("image_map")

    image_numbers = []
    for acquisition_id in acquisition_ids:
        image_number = image_map.get(str(acquisition_id))
        image_numbers.append(image_number)

    if not cell_input or not channel_map or len(image_numbers) == 0:
        raise HTTPException(status_code=400, detail="The dataset does not have a proper input.")

    df = pd.read_feather(cell_input.get("location"))
    df = df[df["ImageNumber"].isin(image_numbers)]

    content = {
        "x": {
            "label": marker_x,
            "data": (df[f"Intensity_MeanIntensity_FullStack_c{channel_map[marker_x]}"] * 2 ** 16).tolist(),
        },
        "y": {
            "label": marker_y,
            "data": (df[f"Intensity_MeanIntensity_FullStack_c{channel_map[marker_y]}"] * 2 ** 16).tolist(),
        },
    }

    if marker_z:
        content["z"] = {
            "label": marker_z,
            "data": (df[f"Intensity_MeanIntensity_FullStack_c{channel_map[marker_z]}"] * 2 ** 16).tolist(),
        }

    if heatmap_type and heatmap:
        if heatmap_type == "channel":
            channel_map = dataset.input.get("channel_map")
            heatmap_data = (df[f"Intensity_MeanIntensity_FullStack_c{channel_map[heatmap]}"] * 2 ** 16).tolist()
        else:
            heatmap_data = df[heatmap].tolist()

        content["heatmap"] = {"label": heatmap, "data": heatmap_data}
    elif len(acquisition_ids) > 1:
        image_map_inv = {v: k for k, v in image_map.items()}
        content["heatmap"] = {
            "label": "Acquisition",
            "data": ([image_map_inv.get(item) for item in df["ImageNumber"].tolist()]),
        }

    # await redis_manager.cache.set(request.url.path, ujson.dumps(content))
    return ORJSONResponse(content)


@router.get("/boxplot", response_model=List[PlotSeriesDto])
async def read_box_plot_data(
    dataset_id: int,
    acquisition_id: int,
    markers: List[str] = Query(None),
    current_user: UserModel = Depends(get_current_active_user),
    db: Session = Depends(get_db),
):
    """
    Read box plot data from the dataset
    """
    # content = await redis_manager.cache.get(request.url.path)
    # if content:
    #     return UJSONResponse(content=ujson.loads(content))

    dataset = dataset_crud.get(db, id=dataset_id)
    cell_input = dataset.input.get("cell")
    channel_map = dataset.input.get("channel_map")
    image_map = dataset.input.get("image_map")
    image_number = image_map.get(str(acquisition_id))
    if not cell_input or not image_number or not channel_map:
        raise HTTPException(status_code=400, detail="The dataset does not have a proper input.")

    content = []
    df = pd.read_feather(cell_input.get("location"))
    df = df[df["ImageNumber"] == image_number]

    for marker in markers:
        content.append(
            {
                "label": marker,
                "data": (df[f"Intensity_MeanIntensity_FullStack_c{channel_map[marker]}"] * 2 ** 16).tolist(),
            }
        )

    # await redis_manager.cache.set(request.url.path, ujson.dumps(content))
    return ORJSONResponse(content)


@router.get("/pca", response_model=PcaDto)
async def read_pca_data(
    dataset_id: int,
    n_components: int,
    acquisition_ids: List[int] = Query(None),
    heatmap_type: Optional[str] = None,
    heatmap: Optional[str] = None,
    markers: List[str] = Query(None),
    current_user: UserModel = Depends(get_current_active_user),
    db: Session = Depends(get_db),
):
    """
    Calculate Principal Component Analysis data for the dataset
    """

    content = pca.process_pca(db, dataset_id, acquisition_ids, n_components, markers, heatmap_type, heatmap)
    return ORJSONResponse(content)


@router.post("/tsne")
def submit_tsne(
    params: TsneSubmissionDto,
    current_user: UserModel = Depends(get_current_active_user),
    db: Session = Depends(get_db),
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


@router.get("/tsne", response_model=TsneDto)
async def read_tsne_data(
    dataset_id: int,
    name: str,
    heatmap_type: Optional[str] = None,
    heatmap: Optional[str] = None,
    current_user: UserModel = Depends(get_current_active_user),
    db: Session = Depends(get_db),
):
    """
    Read t-SNE result data
    """

    content = tsne.get_tsne_result(db, dataset_id, name, heatmap_type, heatmap)
    return ORJSONResponse(content)


@router.post("/umap")
def submit_umap(
    params: UmapSubmissionDto,
    current_user: UserModel = Depends(get_current_active_user),
    db: Session = Depends(get_db),
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


@router.get("/umap", response_model=UmapDto)
async def read_umap_data(
    dataset_id: int,
    name: str,
    heatmap_type: Optional[str] = None,
    heatmap: Optional[str] = None,
    current_user: UserModel = Depends(get_current_active_user),
    db: Session = Depends(get_db),
):
    """
    Read UMAP result data
    """

    content = umap.get_umap_result(db, dataset_id, name, heatmap_type, heatmap)
    return ORJSONResponse(content)


@router.post("/phenograph")
def submit_phenograph(
    params: PhenographSubmissionDto,
    current_user: UserModel = Depends(get_current_active_user),
    db: Session = Depends(get_db),
):
    """
    Start PhenoGraph data processing
    """

    worker.process_phenograph.send(
        params.dataset_id,
        params.acquisition_ids,
        params.markers,
        params.nearest_neighbors,
        params.jaccard,
        params.primary_metric,
        params.min_cluster_size,
    )
    return ORJSONResponse({"status": "submitted"})


@router.get("/phenograph", response_model=PhenographDto)
async def read_phenograph_data(
    dataset_id: int,
    name: str,
    current_user: UserModel = Depends(get_current_active_user),
    db: Session = Depends(get_db),
):
    """
    Read PhenoGraph result data
    """

    content = phenograph.get_phenograph_result(db, dataset_id, name)
    return ORJSONResponse(content)


@router.post("/region/stats", response_model=List[RegionChannelStatsDto])
async def calculate_region_stats(
    params: RegionStatsSubmissionDto,
    request: Request,
    current_user: UserModel = Depends(get_current_active_user),
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
