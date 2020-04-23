import logging
import os
import uuid
from typing import List, Set

from fastapi import APIRouter, Depends, File, HTTPException, UploadFile
from sqlalchemy.orm import Session

import histocat.worker as worker
from histocat.api.utils.db import get_db
from histocat.api.utils.security import get_current_active_user
from histocat.config import config
from histocat.modules.user.models import UserModel

from . import service
from .dto import (
    ExperimentCreateDto,
    ExperimentDatasetDto,
    ExperimentDto,
    ExperimentUpdateDto,
)

logger = logging.getLogger(__name__)
router = APIRouter()


@router.get("/", response_model=List[ExperimentDto])
def read_all(
    db: Session = Depends(get_db), current_user: UserModel = Depends(get_current_active_user),
):
    """
    Retrieve experiments
    """
    items = service.get_multi(db, user=current_user)
    return items


@router.get("/tags", response_model=Set[str])
def read_tags(db: Session = Depends(get_db), current_user: UserModel = Depends(get_current_active_user)):
    """
    Retrieve tags
    """
    items = service.get_tags(db)
    return items


@router.post("/", response_model=ExperimentDto)
def create(
    *,
    db: Session = Depends(get_db),
    params: ExperimentCreateDto,
    current_user: UserModel = Depends(get_current_active_user),
):
    """
    Create new experiment
    """
    if not current_user.is_active:
        raise HTTPException(status_code=401, detail="The user cannot create experiments.")

    item = service.get_by_name(db, name=params.name)
    if item:
        raise HTTPException(
            status_code=400, detail="The experiment with this name already exists in the system.",
        )
    item = service.create(db, user_id=current_user.id, params=params)
    return item


@router.get("/{id}", response_model=ExperimentDto)
def read_by_id(
    id: int, current_user: UserModel = Depends(get_current_active_user), db: Session = Depends(get_db),
):
    """
    Get a specific experiment by id
    """
    item = service.get(db, id=id)
    return item


@router.delete("/{id}", response_model=ExperimentDto)
def delete_by_id(
    id: int, current_user: UserModel = Depends(get_current_active_user), db: Session = Depends(get_db),
):
    """
    Delete a specific experiment by id
    """
    item = service.remove(db, id=id)
    return item


@router.put("/{id}", response_model=ExperimentDto)
def update(
    *,
    db: Session = Depends(get_db),
    id: int,
    params: ExperimentUpdateDto,
    current_user: UserModel = Depends(get_current_active_user),
):
    """
    Update an experiment
    """
    item = service.get(db, id=id)

    if not item:
        raise HTTPException(
            status_code=404, detail="The experiment with this id does not exist in the system",
        )
    item = service.update(db, item=item, params=params)
    return item


@router.post("/{id}/upload")
def upload_data(
    id: int,
    file: UploadFile = File(None),
    user: UserModel = Depends(get_current_active_user),
    db: Session = Depends(get_db),
):
    path = os.path.join(config.INBOX_DIRECTORY, str(uuid.uuid4()))
    if not os.path.exists(path):
        os.makedirs(path)
    uri = os.path.join(path, file.filename)
    with open(uri, "wb") as f:
        f.write(file.file.read())
    worker.import_data.send(uri, id, user.id)
    return {"uri": uri}


@router.get("/{id}/data", response_model=ExperimentDatasetDto)
async def read_data(
    id: int, current_user: UserModel = Depends(get_current_active_user), db: Session = Depends(get_db),
):
    """
    Get all experiment data
    """
    item = service.get_data(db, id=id)
    return item
