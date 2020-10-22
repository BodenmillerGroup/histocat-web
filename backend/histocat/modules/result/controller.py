import logging
import os
from io import BytesIO
from typing import Sequence
from zipfile import ZipFile, ZIP_DEFLATED

import anndata as ad
from fastapi import APIRouter, Depends, HTTPException
from fastapi.responses import ORJSONResponse
from sqlalchemy.orm import Session
from starlette.responses import StreamingResponse

from histocat.api.db import get_db
from histocat.api.security import get_active_member, get_active_user
from histocat.modules.member.models import MemberModel
from histocat.modules.user.models import UserModel

from . import service
from .dto import ResultDto, ResultUpdateDto
from ..analysis.dto import PcaDto, TsneDto, UmapDto
from ...core.utils import stream_bytes
from ...io.dataset import ANNDATA_FILE_EXTENSION

logger = logging.getLogger(__name__)
router = APIRouter()


@router.get("/groups/{group_id}/datasets/{dataset_id}/results", response_model=Sequence[ResultDto])
def get_dataset_results(
    group_id: int, dataset_id: int, db: Session = Depends(get_db), member: MemberModel = Depends(get_active_member),
):
    """
    Retrieve own results for specified dataset
    """
    items = service.get_dataset_results(db, dataset_id=dataset_id)
    return items


@router.get("/groups/{group_id}/results/{result_id}", response_model=ResultDto)
def get_by_id(
    group_id: int, result_id: int, user: UserModel = Depends(get_active_user), db: Session = Depends(get_db),
):
    """
    Get a specific result by id
    """
    item = service.get(db, id=result_id)
    return item


@router.patch("/groups/{group_id}/results/{result_id}", response_model=ResultDto)
def update(
    group_id: int,
    result_id: int,
    params: ResultUpdateDto,
    member: MemberModel = Depends(get_active_member),
    db: Session = Depends(get_db),
):
    """
    Update result
    """
    item = service.get(db, id=result_id)
    if not item:
        raise HTTPException(
            status_code=404, detail="Result not found",
        )
    item = service.update(db, item=item, params=params)
    return item


@router.delete("/groups/{group_id}/results/{result_id}", response_model=ResultDto)
def delete_by_id(
    group_id: int, result_id: int, user: UserModel = Depends(get_active_user), db: Session = Depends(get_db),
):
    """
    Delete a specific dataset by id
    """
    item = service.remove(db, id=result_id)
    return item


@router.get("/results/{result_id}/download")
async def download_by_id(result_id: int, db: Session = Depends(get_db)):
    """
    Download result by id
    """
    item = service.get(db, id=result_id)
    if item is None:
        raise HTTPException(status_code=404, detail=f"Cannot find result [{result_id}]")

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


@router.get("/results/{result_id}/pca", response_model=PcaDto)
async def get_pca_data(
    result_id: int,
    user: UserModel = Depends(get_active_user),
    db: Session = Depends(get_db),
):
    """Get PCA data"""

    result = service.get(db, id=result_id)
    if not result:
        raise HTTPException(status_code=404, detail="Result not found.")
    location = os.path.join(result.location, f"output{ANNDATA_FILE_EXTENSION}")

    adata = ad.read_h5ad(location)

    output = {
        "acquisitionIds": adata.obs["AcquisitionId"].tolist(),
        "cellIds": adata.obs["CellId"].tolist(),
        "objectNumbers": adata.obs["ObjectNumber"].tolist(),
        "x": {"label": "PCA1", "data": adata.obsm["X_pca"][:, 0].tolist(), },
        "y": {"label": "PCA2", "data": adata.obsm["X_pca"][:, 1].tolist(), },
    }

    return ORJSONResponse(output)


@router.get("/results/{result_id}/tsne", response_model=TsneDto)
async def get_tsne_data(
    result_id: int,
    user: UserModel = Depends(get_active_user),
    db: Session = Depends(get_db),
):
    """Get tSNE data"""

    result = service.get(db, id=result_id)
    if not result:
        raise HTTPException(status_code=404, detail="Result not found.")
    location = os.path.join(result.location, f"output{ANNDATA_FILE_EXTENSION}")

    adata = ad.read_h5ad(location)

    output = {
        "acquisitionIds": adata.obs["AcquisitionId"].tolist(),
        "cellIds": adata.obs["CellId"].tolist(),
        "objectNumbers": adata.obs["ObjectNumber"].tolist(),
        "x": {"label": "tSNE1", "data": adata.obsm["X_tsne"][:, 0].tolist(), },
        "y": {"label": "tSNE2", "data": adata.obsm["X_tsne"][:, 1].tolist(), },
    }

    return ORJSONResponse(output)


@router.get("/results/{result_id}/umap", response_model=UmapDto)
async def get_umap_data(
    result_id: int,
    user: UserModel = Depends(get_active_user),
    db: Session = Depends(get_db),
):
    """Get UMAP data"""

    result = service.get(db, id=result_id)
    if not result:
        raise HTTPException(status_code=404, detail="Result not found.")
    location = os.path.join(result.location, f"output{ANNDATA_FILE_EXTENSION}")

    adata = ad.read_h5ad(location)

    output = {
        "acquisitionIds": adata.obs["AcquisitionId"].tolist(),
        "cellIds": adata.obs["CellId"].tolist(),
        "objectNumbers": adata.obs["ObjectNumber"].tolist(),
        "x": {"label": "tSNE1", "data": adata.obsm["X_umap"][:, 0].tolist(), },
        "y": {"label": "tSNE2", "data": adata.obsm["X_umap"][:, 1].tolist(), },
    }

    return ORJSONResponse(output)
