import logging
import os
from io import BytesIO
from typing import Optional, Sequence
from zipfile import ZIP_DEFLATED, ZipFile

import anndata as ad
import matplotlib
import scanpy as sc
from fastapi import APIRouter, Depends, HTTPException
from fastapi.responses import FileResponse, ORJSONResponse, StreamingResponse
from sqlalchemy.orm import Session
from starlette.status import HTTP_404_NOT_FOUND

from histocat.api.db import get_db
from histocat.api.security import get_active_member, get_active_user
from histocat.core.constants import ANNDATA_FILE_EXTENSION
from histocat.core.image import get_qualitative_colors, get_sequential_colors
from histocat.core.member.models import MemberModel
from histocat.core.result import service
from histocat.core.result.dto import (
    ColorsDataDto,
    ResultDataDto,
    ResultDto,
    ResultUpdateDto,
)
from histocat.core.user.models import UserModel
from histocat.core.utils import stream_bytes

logger = logging.getLogger(__name__)
router = APIRouter()


@router.get("/groups/{group_id}/datasets/{dataset_id}/results", response_model=Sequence[ResultDto])
def get_dataset_results(
    group_id: int,
    dataset_id: int,
    db: Session = Depends(get_db),
    member: MemberModel = Depends(get_active_member),
):
    """Retrieve own results for specified dataset"""
    items = service.get_dataset_results(db, dataset_id=dataset_id)
    return items


@router.get("/groups/{group_id}/results/{result_id}", response_model=ResultDto)
def get_by_id(
    group_id: int,
    result_id: int,
    member: MemberModel = Depends(get_active_member),
    db: Session = Depends(get_db),
):
    """Get result by id"""
    item = service.get(db, id=result_id)
    if not item:
        raise HTTPException(status_code=HTTP_404_NOT_FOUND, detail=f"Result id:{result_id} not found")
    return item


@router.get("/groups/{group_id}/results/{result_id}/data", response_model=ResultDataDto)
def get_result_data(
    group_id: int,
    result_id: int,
    member: MemberModel = Depends(get_active_member),
    db: Session = Depends(get_db),
):
    """Get result data by id"""
    result = service.get(db, id=result_id)
    if not result:
        raise HTTPException(status_code=HTTP_404_NOT_FOUND, detail=f"Result id:{result_id} not found")

    location = os.path.join(result.location, f"output{ANNDATA_FILE_EXTENSION}")
    # Read AnnData file
    adata = ad.read_h5ad(location)

    output = {
        "cellIds": adata.obs["CellId"].tolist(),
        "markers": adata.var_names.tolist(),
    }

    mappings = {}
    if "pca" in result.output:
        mappings["pca"] = {
            "x": {
                "label": "PCA1",
                "data": adata.obsm["X_pca"][:, 0].tolist(),
            },
            "y": {
                "label": "PCA2",
                "data": adata.obsm["X_pca"][:, 1].tolist(),
            },
        }

    if "tsne" in result.output:
        mappings["tsne"] = {
            "x": {
                "label": "tSNE1",
                "data": adata.obsm["X_tsne"][:, 0].tolist(),
            },
            "y": {
                "label": "tSNE2",
                "data": adata.obsm["X_tsne"][:, 1].tolist(),
            },
        }

    if "umap" in result.output:
        mappings["umap"] = {
            "x": {
                "label": "UMAP1",
                "data": adata.obsm["X_umap"][:, 0].tolist(),
            },
            "y": {
                "label": "UMAP2",
                "data": adata.obsm["X_umap"][:, 1].tolist(),
            },
        }

    output["mappings"] = mappings
    return ORJSONResponse(output)


@router.get("/groups/{group_id}/results/{result_id}/colors", response_model=ColorsDataDto)
def get_colors_data(
    group_id: int,
    result_id: int,
    colors_type: str,
    colors_name: str,
    member: MemberModel = Depends(get_active_member),
    db: Session = Depends(get_db),
):
    """Get colors (heatmap) data"""
    result = service.get(db, id=result_id)
    if not result:
        raise HTTPException(status_code=HTTP_404_NOT_FOUND, detail=f"Result id:{result_id} not found")

    location = os.path.join(result.location, f"output{ANNDATA_FILE_EXTENSION}")
    # Read AnnData file
    adata = ad.read_h5ad(location)

    output = {
        "cellIds": adata.obs["CellId"].tolist(),
    }

    colors = None
    if colors_type == "marker":
        values = adata.X[:, adata.var.index == colors_name]
        mappable = get_sequential_colors()
        colors = [matplotlib.colors.rgb2hex(c) for c in mappable.to_rgba(values)]
    elif colors_type == "clustering":
        values = sc.get.obs_df(adata, keys=[colors_name]).astype(int)
        mappable = get_qualitative_colors()
        colors = [matplotlib.colors.rgb2hex(c) for c in mappable.to_rgba(values)]

    if colors is not None:
        output["colors"] = {"type": colors_type, "name": colors_name, "data": colors}

    return ORJSONResponse(output)


@router.patch("/groups/{group_id}/results/{result_id}", response_model=ResultDto)
def update(
    group_id: int,
    result_id: int,
    params: ResultUpdateDto,
    member: MemberModel = Depends(get_active_member),
    db: Session = Depends(get_db),
):
    """Update result"""
    item = service.get(db, id=result_id)
    if not item:
        raise HTTPException(status_code=HTTP_404_NOT_FOUND, detail=f"Result id:{result_id} not found")
    item = service.update(db, item=item, params=params)
    return item


@router.delete("/groups/{group_id}/results/{result_id}", response_model=ResultDto)
def delete_by_id(
    group_id: int,
    result_id: int,
    user: UserModel = Depends(get_active_user),
    db: Session = Depends(get_db),
):
    """Delete a specific dataset by id"""
    item = service.remove(db, id=result_id)
    return item


@router.get("/results/{result_id}/download")
async def download_by_id(result_id: int, db: Session = Depends(get_db)):
    """Download result by id"""
    item = service.get(db, id=result_id)
    if not item:
        raise HTTPException(status_code=HTTP_404_NOT_FOUND, detail=f"Result id:{result_id} not found")

    file_name = f"{item.name}.zip"
    abs_src = os.path.abspath(item.location)
    buffer = BytesIO()
    with ZipFile(buffer, "w", ZIP_DEFLATED) as zip:
        for folderName, _, filenames in os.walk(item.location):
            for filename in filenames:
                absname = os.path.abspath(os.path.join(folderName, filename))
                arcname = absname[len(abs_src) + 1 :]
                zip.write(absname, arcname)

    headers = {"Content-Disposition": f'attachment; filename="{file_name}"'}
    return StreamingResponse(stream_bytes(buffer.getvalue()), media_type="application/zip", headers=headers)


@router.get("/results/{result_id}/plot", responses={200: {"content": {"image/png": {}}}})
async def get_plot_image(
    result_id: int,
    plot_name: Optional[str] = None,
    # user: UserModel = Depends(get_current_active_user),
    db: Session = Depends(get_db),
):
    """Get plot image by name"""
    result = service.get(db, id=result_id)
    if not result:
        raise HTTPException(status_code=HTTP_404_NOT_FOUND, detail=f"Result id:{result_id} not found")

    filename = f"{plot_name}.png"
    location = os.path.join(result.location, filename)

    return FileResponse(location, media_type="image/png", filename=filename)


@router.get("/groups/{group_id}/results/{result_id}/scatterplot")
async def get_scatter_plot_data(
    group_id: int,
    result_id: int,
    marker_x: str,
    marker_y: str,
    member: MemberModel = Depends(get_active_member),
    db: Session = Depends(get_db),
):
    """Read scatter plot data from the result."""

    result = service.get(db, id=result_id)
    if not result:
        raise HTTPException(status_code=HTTP_404_NOT_FOUND, detail=f"Result id:{result_id} not found")
    location = os.path.join(result.location, f"output{ANNDATA_FILE_EXTENSION}")

    adata = ad.read_h5ad(location)

    output = {
        "cellIds": adata.obs["CellId"].tolist(),
        "x": {
            "label": marker_x,
            "data": adata.X[:, adata.var.index == marker_x][:, 0].tolist(),
        },
        "y": {
            "label": marker_y,
            "data": adata.X[:, adata.var.index == marker_y][:, 0].tolist(),
        },
    }

    return ORJSONResponse(output)
