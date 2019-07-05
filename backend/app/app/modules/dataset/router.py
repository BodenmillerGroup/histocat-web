from typing import List

from fastapi import APIRouter, Depends
from sqlalchemy.orm import Session

from app.api.utils.db import get_db
from app.api.utils.security import get_current_active_superuser, get_current_active_user
from app.modules.user.db import User
from . import crud
from .models import DatasetModel, DatasetCreateModel
import app.worker as worker

router = APIRouter()


@router.get("/", response_model=List[DatasetModel])
def read_all(
        db: Session = Depends(get_db),
        skip: int = 0,
        limit: int = 100,
        current_user: User = Depends(get_current_active_superuser),
):
    """
    Retrieve datasets
    """
    items = crud.get_multi(db, skip=skip, limit=limit)
    return items


@router.get("/experiment/{experiment_id}", response_model=List[DatasetModel])
def read_own_by_experiment(
        experiment_id: int,
        db: Session = Depends(get_db),
        current_user: User = Depends(get_current_active_superuser),
):
    """
    Retrieve own datasets for specified experiment
    """
    items = crud.get_own_by_experiment_id(db, user_id=current_user.id, experiment_id=experiment_id)
    return items


@router.get("/{id}", response_model=DatasetModel)
def read_by_id(
        id: int,
        current_user: User = Depends(get_current_active_user),
        db: Session = Depends(get_db),
):
    """
    Get a specific dataset by id
    """
    item = crud.get(db, id=id)
    return item


@router.post("/", response_model=DatasetModel)
def create(
        *,
        db: Session = Depends(get_db),
        params: DatasetCreateModel,
        current_user: User = Depends(get_current_active_superuser),
):
    """
    Create new dataset
    """
    item = crud.create(db, user_id=current_user.id, params=params)
    worker.prepare_dataset.send(item.id)
    return item


@router.delete("/{id}", response_model=DatasetModel)
def delete_by_id(
        id: int,
        current_user: User = Depends(get_current_active_superuser),
        db: Session = Depends(get_db),
):
    """
    Delete a specific dataset by id
    """
    item = crud.remove(db, id=id)
    return item
