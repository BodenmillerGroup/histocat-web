from typing import List

from fastapi import APIRouter, Depends, HTTPException
from sqlalchemy.orm import Session

from app import crud
from app.api.utils.db import get_db
from app.api.utils.security import get_current_active_superuser, get_current_active_user
from app.db_models import User
from app.models.experiment import ExperimentModel, ExperimentInCreateModel, ExperimentInUpdateModel

router = APIRouter()


@router.get("/experiments/", tags=["experiments"], response_model=List[ExperimentModel])
def read_experiments(
    db: Session = Depends(get_db),
    skip: int = 0,
    limit: int = 100,
    current_user: User = Depends(get_current_active_superuser),
):
    """
    Retrieve experiments
    """
    expperiments = crud.experiment.get_multi(db, skip=skip, limit=limit)
    return expperiments


@router.post("/experiments/", tags=["experiments"], response_model=ExperimentModel)
def create_experiment(
    *,
    db: Session = Depends(get_db),
    experiment_in: ExperimentInCreateModel,
    current_user: User = Depends(get_current_active_superuser),
):
    """
    Create new experiment
    """
    experiment = crud.experiment.get_by_name(db, name=experiment_in.name)
    if experiment:
        raise HTTPException(
            status_code=400,
            detail="The experiment with this name already exists in the system.",
        )
    experiment = crud.experiment.create(db, params=experiment_in)
    return experiment


@router.get("/experiments/{experiment_id}", tags=["experiments"], response_model=ExperimentModel)
def read_experiment_by_id(
    experiment_id: int,
    current_user: User = Depends(get_current_active_user),
    db: Session = Depends(get_db),
):
    """
    Get a specific experiment by id
    """
    experiment = crud.experiment.get(db, id=experiment_id)
    return experiment


@router.put("/experiments/{experiment_id}", tags=["experiments"], response_model=ExperimentModel)
def update_experiment(
    *,
    db: Session = Depends(get_db),
    experiment_id: int,
    experiment_in: ExperimentInUpdateModel,
    current_user: User = Depends(get_current_active_superuser),
):
    """
    Update a user
    """
    experiment = crud.experiment.get(db, id=experiment_id)

    if not experiment:
        raise HTTPException(
            status_code=404,
            detail="The experiment with this id does not exist in the system",
        )
    experiment = crud.experiment.update(db, experiment=experiment, params=experiment_in)
    return experiment
