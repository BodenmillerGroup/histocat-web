import logging
import os
from io import BytesIO
from typing import List
from zipfile import ZIP_DEFLATED, ZipFile

from fastapi import APIRouter, Depends, HTTPException
from fastapi.responses import ORJSONResponse
from sqlalchemy.orm import Session
from starlette.responses import StreamingResponse

from app.api.utils.db import get_db
from app.api.utils.security import get_current_active_user
from app.core.utils import stream_bytes
from app.modules.user.models import UserModel

from . import service
from .dto import DatasetDto

logger = logging.getLogger(__name__)
router = APIRouter()


@router.get("/", response_model=List[DatasetDto], response_class=ORJSONResponse)
def read_all(
    db: Session = Depends(get_db), current_user: UserModel = Depends(get_current_active_user),
):
    """
    Retrieve datasets
    """
    items = service.get_multi(db)
    return items


@router.get("/experiment/{experiment_id}", response_model=List[DatasetDto], response_class=ORJSONResponse)
def read_own_by_experiment(
    experiment_id: int, db: Session = Depends(get_db), current_user: UserModel = Depends(get_current_active_user),
):
    """
    Retrieve own datasets for specified experiment
    """
    items = service.get_own_by_experiment_id(db, experiment_id=experiment_id)
    return items


@router.get("/{id}", response_model=DatasetDto, response_class=ORJSONResponse)
def read_by_id(
    id: int, current_user: UserModel = Depends(get_current_active_user), db: Session = Depends(get_db),
):
    """
    Get a specific dataset by id
    """
    item = service.get(db, id=id)
    return item


@router.delete("/{id}", response_model=DatasetDto, response_class=ORJSONResponse)
def delete_by_id(
    id: int, current_user: UserModel = Depends(get_current_active_user), db: Session = Depends(get_db),
):
    """
    Delete a specific dataset by id
    """
    item = service.remove(db, id=id)
    return item


@router.get("/{id}/download")
async def download_by_id(id: int, db: Session = Depends(get_db)):
    """
    Download dataset by id
    """
    item = service.get(db, id=id)
    if item is None:
        raise HTTPException(status_code=400, detail=f"Cannot find dataset [{id}]")

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
