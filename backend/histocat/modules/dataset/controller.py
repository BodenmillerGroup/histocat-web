import logging
import os
from io import BytesIO
from typing import Sequence
from zipfile import ZIP_DEFLATED, ZipFile

from fastapi import APIRouter, Depends, HTTPException
from fastapi.responses import ORJSONResponse
from sqlalchemy.orm import Session
from starlette.responses import StreamingResponse

from histocat.api.db import get_db
from histocat.api.security import get_active_member, get_active_user
from histocat.core.utils import stream_bytes
from histocat.modules.member.models import MemberModel
from histocat.modules.user.models import UserModel

from . import service
from .dto import DatasetDto, DatasetUpdateDto

logger = logging.getLogger(__name__)
router = APIRouter()


@router.get("/groups/{group_id}/projects/{project_id}/datasets", response_model=Sequence[DatasetDto])
def get_project_datasets(
    group_id: int, project_id: int, db: Session = Depends(get_db), member: MemberModel = Depends(get_active_member),
):
    """
    Retrieve own project's datasets
    """
    items = service.get_project_datasets(db, project_id=project_id)
    return items


@router.patch("/groups/{group_id}/datasets/{dataset_id}", response_model=DatasetDto)
def update(
    group_id: int,
    dataset_id: int,
    params: DatasetUpdateDto,
    member: MemberModel = Depends(get_active_member),
    db: Session = Depends(get_db),
):
    """
    Update dataset
    """
    item = service.get(db, id=dataset_id)
    if not item:
        raise HTTPException(
            status_code=404, detail="Dataset not found",
        )
    item = service.update(db, item=item, params=params)
    return item


@router.get("/datasets/{dataset_id}/centroids")
def get_centroids(
    dataset_id: int, user: UserModel = Depends(get_active_user), db: Session = Depends(get_db),
):
    """
    Get dataset cell centroids
    """
    dataset = service.get(db, id=dataset_id)
    if dataset is None:
        raise HTTPException(status_code=400, detail=f"Cannot find dataset [{dataset_id}]")
    content = service.get_centroids(dataset)
    return ORJSONResponse(content)


@router.get("/datasets/{dataset_id}", response_model=DatasetDto)
def get_by_id(
    dataset_id: int, user: UserModel = Depends(get_active_user), db: Session = Depends(get_db),
):
    """
    Get a specific dataset by id
    """
    item = service.get(db, id=dataset_id)
    return item


@router.delete("/datasets/{dataset_id}", response_model=DatasetDto)
def delete_by_id(
    dataset_id: int, user: UserModel = Depends(get_active_user), db: Session = Depends(get_db),
):
    """
    Delete a specific dataset by id
    """
    item = service.remove(db, id=dataset_id)
    return item


@router.get("/datasets/{dataset_id}/download")
async def download_by_id(dataset_id: int, db: Session = Depends(get_db)):
    """
    Download dataset by id
    """
    item = service.get(db, id=dataset_id)
    if item is None:
        raise HTTPException(status_code=400, detail=f"Cannot find dataset [{dataset_id}]")

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
