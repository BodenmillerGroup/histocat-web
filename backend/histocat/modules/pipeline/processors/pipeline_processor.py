from typing import Sequence, Any
import os
from datetime import datetime

import scanpy as sc
from fastapi import HTTPException
from fastapi.encoders import jsonable_encoder
from sqlalchemy.orm import Session

from histocat.core.utils import timeit
from histocat.core.notifier import Message
from histocat.core.redis_manager import UPDATES_CHANNEL_NAME, redis_manager
from histocat.io.dataset import ANNDATA_FILE_EXTENSION
from histocat.modules.result import service as result_service
from histocat.modules.dataset import service as dataset_service
from histocat.modules.pipeline.processors.steps.filter_acquisitions import filter_acquisitions
from histocat.modules.result.dto import ResultCreateDto


@timeit
def process_pipeline(
    db: Session,
    dataset_id: int,
    acquisition_ids: Sequence[int],
    steps: Sequence[Any]
):
    """
    Process pipeline steps
    """

    dataset = dataset_service.get(db, id=dataset_id)
    cell_input = dataset.meta.get("cell")

    if not cell_input or len(acquisition_ids) == 0:
        raise HTTPException(status_code=400, detail="The dataset does not have a proper input.")

    adata = sc.read_h5ad(cell_input.get("location"))

    # Subset observations for selected acquisitions
    adata = filter_acquisitions(adata, acquisition_ids)

    result_create_params = ResultCreateDto(
        dataset_id=dataset_id,
        status="ready",
        name=str(datetime.utcnow()),
        input=acquisition_ids,
        pipeline=steps
    )
    result = result_service.create(db, params=result_create_params)

    location = os.path.join(result.location, f"output{ANNDATA_FILE_EXTENSION}")
    sc.write(location, adata=adata)

    redis_manager.publish(
        UPDATES_CHANNEL_NAME, Message(dataset.project_id, "result_ready", jsonable_encoder(result)),
    )
