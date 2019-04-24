import logging
from typing import List

from fastapi import APIRouter, Depends, HTTPException, UploadFile, File
from sqlalchemy.orm import Session
from starlette.background import BackgroundTasks

from app.api.utils.db import get_db
from app.api.utils.security import get_current_active_superuser, get_current_active_user
from app.io.mcd_loader import McdLoader
from app.modules.user.db import User
from . import crud
from .models import ExperimentModel, ExperimentCreateModel, ExperimentUpdateModel

logger = logging.getLogger(__name__)

router = APIRouter()


@router.get("/", response_model=List[ExperimentModel])
def read_experiments(
    db: Session = Depends(get_db),
    skip: int = 0,
    limit: int = 100,
    current_user: User = Depends(get_current_active_superuser),
):
    """
    Retrieve experiments
    """
    items = crud.get_multi(db, skip=skip, limit=limit)
    return items


@router.post("/", response_model=ExperimentModel)
def create_experiment(
    *,
    db: Session = Depends(get_db),
    params: ExperimentCreateModel,
    current_user: User = Depends(get_current_active_superuser),
):
    """
    Create new experiment
    """
    item = crud.get_by_name(db, name=params.name)
    if item:
        raise HTTPException(
            status_code=400,
            detail="The experiment with this name already exists in the system.",
        )
    item = crud.create(db, params=params)
    return item


@router.get("/{id}", response_model=ExperimentModel)
def read_experiment_by_id(
    id: int,
    current_user: User = Depends(get_current_active_user),
    db: Session = Depends(get_db),
):
    """
    Get a specific experiment by id
    """
    item = crud.get(db, id=id)
    return item


@router.put("/{id}", response_model=ExperimentModel)
def update_experiment(
    *,
    db: Session = Depends(get_db),
    id: int,
    params: ExperimentUpdateModel,
    current_user: User = Depends(get_current_active_superuser),
):
    """
    Update an experiment
    """
    item = crud.get(db, id=id)

    if not item:
        raise HTTPException(
            status_code=404,
            detail="The experiment with this id does not exist in the system",
        )
    item = crud.update(db, item=item, params=params)
    return item


@router.post("/{id}/upload_slide")
def upload_slide(id: int, background_tasks: BackgroundTasks, file: UploadFile = File(None), db: Session = Depends(get_db)):
    experiment = crud.get(db, id=id)
    background_tasks.add_task(McdLoader.load, file, db, experiment)
    # await McdLoader.load(file, db, experiment)
    return {"filename": file.filename}
