from typing import List, Optional

from fastapi.encoders import jsonable_encoder

from app.modules.experiment.db import Experiment
from .models import ExperimentInCreateModel, ExperimentInUpdateModel


def get(db_session, *, id: int) -> Optional[Experiment]:
    return db_session.query(Experiment).filter(Experiment.id == id).first()


def get_by_name(db_session, *, name: str) -> Optional[Experiment]:
    return db_session.query(Experiment).filter(Experiment.name == name).first()


def get_multi(db_session, *, skip: int = 0, limit: int = 100) -> List[Optional[Experiment]]:
    return db_session.query(Experiment).offset(skip).limit(limit).all()


def create(db_session, *, params: ExperimentInCreateModel) -> Experiment:
    experiment = Experiment(
        name=params.name,
        root_directory=params.root_directory,
        description=params.description
    )
    db_session.add(experiment)
    db_session.commit()
    db_session.refresh(experiment)
    return experiment


def update(db_session, *, experiment: Experiment, params: ExperimentInUpdateModel) -> Experiment:
    experiment_data = jsonable_encoder(experiment)
    for field in experiment_data:
        if field in params.fields:
            value_in = getattr(params, field)
            if value_in is not None:
                setattr(experiment, field, value_in)
    db_session.add(experiment)
    db_session.commit()
    db_session.refresh(experiment)
    return experiment
