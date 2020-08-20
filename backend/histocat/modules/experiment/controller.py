import logging
import os
import uuid
from typing import Sequence, Set

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


@router.get("/groups/{id}/tags", response_model=Set[str])
def get_tags(id: int, db: Session = Depends(get_db), user: UserModel = Depends(get_current_active_user)):
    """
    Get group experiments tags
    """
    items = service.get_tags(db, group_id=id)
    return items


@router.get("/groups/{id}/experiments", response_model=Sequence[ExperimentDto])
def get_group_experiments(id: int, db: Session = Depends(get_db), user: UserModel = Depends(get_current_active_user)):
    """
    Get group experiments
    """
    items = service.get_group_experiments(db, group_id=id)
    return items


@router.post("/groups/{id}/experiments", response_model=ExperimentDto)
def create(
    id: int,
    params: ExperimentCreateDto,
    db: Session = Depends(get_db),
    user: UserModel = Depends(get_current_active_user),
):
    """
    Create new experiment
    """
    if not user.is_active:
        raise HTTPException(status_code=401, detail="The user cannot create experiments.")

    item = service.get_by_name(db, name=params.name)
    if item:
        raise HTTPException(
            status_code=400, detail="The experiment with this name already exists in the system.",
        )
    item = service.create(db, group_id=id, params=params)
    return item


@router.get("/experiments/{experiment_id}", response_model=ExperimentDto)
def get_by_id(
    experiment_id: int, user: UserModel = Depends(get_current_active_user), db: Session = Depends(get_db),
):
    """
    Get experiment by id
    """
    item = service.get(db, id=experiment_id)
    return item


@router.delete("/experiments/{experiment_id}", response_model=ExperimentDto)
def delete_by_id(
    experiment_id: int, user: UserModel = Depends(get_current_active_user), db: Session = Depends(get_db),
):
    """
    Delete a specific experiment by id
    """
    item = service.remove(db, id=experiment_id)
    return item


@router.put("/experiments/{experiment_id}", response_model=ExperimentDto)
def update(
    *,
    db: Session = Depends(get_db),
    experiment_id: int,
    params: ExperimentUpdateDto,
    user: UserModel = Depends(get_current_active_user),
):
    """
    Update an experiment
    """
    item = service.get(db, id=experiment_id)
    if not item:
        raise HTTPException(
            status_code=404, detail="The experiment with this id does not exist in the system",
        )
    item = service.update(db, item=item, params=params)
    return item


@router.post("/experiments/{experiment_id}/upload")
def upload_data(
    experiment_id: int,
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
    worker.import_data.send(uri, experiment_id, user.id)
    return {"uri": uri}


@router.get("/experiments/{experiment_id}/data", response_model=ExperimentDatasetDto)
async def read_data(
    experiment_id: int, user: UserModel = Depends(get_current_active_user), db: Session = Depends(get_db),
):
    """
    Get all experiment data
    """
    item = service.get_data(db, id=experiment_id)
    return item
