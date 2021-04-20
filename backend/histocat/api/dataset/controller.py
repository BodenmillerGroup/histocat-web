import logging
import os
import uuid
from io import BytesIO
from typing import Sequence
from zipfile import ZIP_DEFLATED, ZipFile

import anndata as ad
import dramatiq
import matplotlib
from fastapi import APIRouter, Depends, File, HTTPException, UploadFile
from fastapi.responses import ORJSONResponse
from sqlalchemy.orm import Session
from starlette.responses import StreamingResponse
from starlette.status import HTTP_404_NOT_FOUND

from histocat.api.db import get_db
from histocat.api.security import get_active_member
from histocat.config import config
from histocat.core.dataset import service as dataset_service
from histocat.core.dataset.dto import DatasetDto, DatasetUpdateDto
from histocat.core.image import get_qualitative_colors
from histocat.core.member.models import MemberModel
from histocat.core.utils import stream_bytes

logger = logging.getLogger(__name__)
router = APIRouter()


@router.get("/groups/{group_id}/projects/{project_id}/datasets", response_model=Sequence[DatasetDto])
def get_project_datasets(
    group_id: int,
    project_id: int,
    db: Session = Depends(get_db),
    member: MemberModel = Depends(get_active_member),
):
    """Retrieve project's datasets"""
    items = dataset_service.get_project_datasets(db, project_id=project_id)
    return items


@router.patch("/groups/{group_id}/datasets/{dataset_id}", response_model=DatasetDto)
def update(
    group_id: int,
    dataset_id: int,
    params: DatasetUpdateDto,
    member: MemberModel = Depends(get_active_member),
    db: Session = Depends(get_db),
):
    """Update dataset"""
    item = dataset_service.get(db, id=dataset_id)
    if not item:
        raise HTTPException(status_code=HTTP_404_NOT_FOUND, detail=f"Dataset id:{dataset_id} not found")
    item = dataset_service.update(db, item=item, params=params)
    return item


@router.get("/groups/{group_id}/datasets/{dataset_id}/centroids")
def get_centroids(
    group_id: int,
    dataset_id: int,
    member: MemberModel = Depends(get_active_member),
    db: Session = Depends(get_db),
):
    """Get dataset cell centroids"""
    dataset = dataset_service.get(db, id=dataset_id)
    if not dataset:
        raise HTTPException(status_code=HTTP_404_NOT_FOUND, detail=f"Dataset id:{dataset_id} not found")

    adata = ad.read_h5ad(dataset.cell_file_location())
    mappable = get_qualitative_colors()
    colors = [matplotlib.colors.rgb2hex(c) for c in mappable.to_rgba(adata.obs["AcquisitionId"])]
    output = {
        "acquisitionIds": adata.obs["AcquisitionId"].tolist(),
        "cellIds": adata.obs["CellId"].tolist(),
        "objectNumbers": adata.obs["ObjectNumber"].tolist(),
        "x": adata.obs["CentroidX"].round(2).tolist(),
        "y": adata.obs["CentroidY"].round(2).tolist(),
        "colors": colors
    }
    return ORJSONResponse(output)


@router.get("/groups/{group_id}/datasets/{dataset_id}", response_model=DatasetDto)
def get_by_id(
    group_id: int,
    dataset_id: int,
    member: MemberModel = Depends(get_active_member),
    db: Session = Depends(get_db),
):
    """Get dataset by id"""
    item = dataset_service.get(db, id=dataset_id)
    if not item:
        raise HTTPException(status_code=HTTP_404_NOT_FOUND, detail=f"Dataset id:{dataset_id} not found")
    return item


@router.delete("/groups/{group_id}/datasets/{dataset_id}", response_model=DatasetDto)
def delete_by_id(
    group_id: int,
    dataset_id: int,
    member: MemberModel = Depends(get_active_member),
    db: Session = Depends(get_db),
):
    """Delete a specific dataset by id"""
    item = dataset_service.remove(db, id=dataset_id)
    if not item:
        raise HTTPException(status_code=HTTP_404_NOT_FOUND, detail=f"Dataset id:{dataset_id} not found")
    return item


@router.get("/datasets/{dataset_id}/download")
async def download_by_id(dataset_id: int, db: Session = Depends(get_db)):
    """Download dataset by id"""
    item = dataset_service.get(db, id=dataset_id)
    if not item:
        raise HTTPException(status_code=HTTP_404_NOT_FOUND, detail=f"Dataset id:{dataset_id} not found")

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


@router.post("/groups/{group_id}/projects/{project_id}/datasets/upload")
def upload_dataset(
    group_id: int,
    project_id: int,
    file: UploadFile = File(None),
    member: MemberModel = Depends(get_active_member),
    db: Session = Depends(get_db),
):
    path = os.path.join(config.INBOX_DIRECTORY, str(uuid.uuid4()))
    if not os.path.exists(path):
        os.makedirs(path)
    uri = os.path.join(path, file.filename)
    with open(uri, "wb") as f:
        f.write(file.file.read())

    broker = dramatiq.get_broker()
    message = dramatiq.Message(
        actor_name="import_dataset",
        queue_name="import",
        args=(),
        kwargs={"uri": uri, "project_id": project_id},
        options={},
    )
    broker.enqueue(message)
    return {"uri": uri}
