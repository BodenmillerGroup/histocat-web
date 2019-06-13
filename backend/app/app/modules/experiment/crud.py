from typing import List, Optional

from fastapi.encoders import jsonable_encoder
from sqlalchemy.orm import Session

from .db import Experiment
from .models import ExperimentCreateModel, ExperimentUpdateModel

ROOT_DIRECTORY = "/data/experiments/"


def get(session: Session, *, id: int) -> Optional[Experiment]:
    return session.query(Experiment).filter(Experiment.id == id).first()


def get_by_name(session: Session, *, name: str) -> Optional[Experiment]:
    return session.query(Experiment).filter(Experiment.name == name).first()


def get_multi(
    session: Session, *, skip: int = 0, limit: int = 100
) -> List[Optional[Experiment]]:
    items = session.query(Experiment).offset(skip).limit(limit).all()
    # TODO: hack to get location property value in API response
    for i in items:
        i.location
    return items


def create(session: Session, *, params: ExperimentCreateModel) -> Experiment:
    entity = Experiment(
        root_directory=ROOT_DIRECTORY,
        name=params.name,
        description=params.description,
        meta=params.meta,
    )
    session.add(entity)
    session.commit()
    session.refresh(entity)
    return entity


def update(
    session: Session, *, item: Experiment, params: ExperimentUpdateModel
) -> Experiment:
    data = jsonable_encoder(item)
    update_data = params.dict(skip_defaults=True)
    for field in data:
        if field in update_data:
            setattr(item, field, update_data[field])
    session.add(item)
    session.commit()
    session.refresh(item)
    return item


def remove(session: Session, *, id: int):
    item = session.query(Experiment).filter(Experiment.id == id).first()
    session.delete(item)
    session.commit()
    return item


def get_dataset(session: Session, *, id: int) -> Optional[Experiment]:
    return session.query(Experiment).filter(Experiment.id == id).first()
