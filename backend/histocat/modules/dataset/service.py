import logging
import os
from typing import List, Optional

from fastapi import HTTPException
from fastapi.encoders import jsonable_encoder
from sqlalchemy.orm import Session
import anndata as ad

from .dto import DatasetCreateDto, DatasetUpdateDto
from .models import DatasetModel

logger = logging.getLogger(__name__)


def get(session: Session, *, id: int) -> Optional[DatasetModel]:
    return session.query(DatasetModel).filter(DatasetModel.id == id).first()


def get_project_datasets(session: Session, *, project_id: int) -> List[DatasetModel]:
    return session.query(DatasetModel).filter(DatasetModel.project_id == project_id).all()


def get_multi(session: Session, *, skip: int = 0, limit: int = 1000) -> List[DatasetModel]:
    return session.query(DatasetModel).offset(skip).limit(limit).all()


def create(session: Session, *, params: DatasetCreateDto) -> DatasetModel:
    data = jsonable_encoder(params)
    entity = DatasetModel(**data)
    session.add(entity)
    session.commit()
    session.refresh(entity)

    entity.location = os.path.join(entity.project.datasets_location, str(entity.id))
    if not os.path.exists(entity.location):
        logger.debug(f"Create location for dataset {entity.id}: {entity.location}")
        os.makedirs(entity.location)

    session.commit()
    session.refresh(entity)
    return entity


def update(session: Session, *, item: DatasetModel, params: DatasetUpdateDto) -> DatasetModel:
    data = item.as_dict()
    update_data = params.dict(exclude_unset=True)
    for field in data:
        if field in update_data:
            setattr(item, field, update_data[field])
    session.add(item)
    session.commit()
    session.refresh(item)
    return item


def remove(session: Session, *, id: int):
    item = session.query(DatasetModel).filter(DatasetModel.id == id).first()
    if item:
        session.delete(item)
        session.commit()
        return item


def get_centroids(dataset: DatasetModel):
    cell_input = dataset.meta.get("cell")

    if not cell_input:
        raise HTTPException(status_code=400, detail="The dataset does not have a proper input.")

    adata = ad.read_h5ad(cell_input.get("location"))
    output = {
        "acquisitionIds": adata.obs["AcquisitionId"].tolist(),
        "cellIds": adata.obs["CellId"].tolist(),
        "objectNumbers": adata.obs["ObjectNumber"].tolist(),
        "x": adata.obs["CentroidX"].round(2).tolist(),
        "y": adata.obs["CentroidY"].round(2).tolist(),
    }
    return output
