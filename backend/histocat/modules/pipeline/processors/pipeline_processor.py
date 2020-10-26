import os
from datetime import datetime
from typing import Any, Sequence

import scanpy as sc
from fastapi import HTTPException
from fastapi.encoders import jsonable_encoder
from sqlalchemy.orm import Session

from histocat.core.notifier import Message
from histocat.core.redis_manager import UPDATES_CHANNEL_NAME, redis_manager
from histocat.core.utils import timeit
from histocat.io.dataset import ANNDATA_FILE_EXTENSION
from histocat.modules.dataset import service as dataset_service
from histocat.modules.pipeline.processors.steps import (
    acquisitions_filter,
    leiden,
    louvain,
    markers_filter,
    neighbors,
    pca,
    scale,
    transformation,
    tsne,
    umap,
)
from histocat.modules.result import service as result_service
from histocat.modules.result.dto import ResultCreateDto, ResultUpdateDto

processors = {
    "markersFilter": markers_filter,
    "transformation": transformation,
    "scale": scale,
    "neighbors": neighbors,
    "pca": pca,
    "tsne": tsne,
    "umap": umap,
    "leiden": leiden,
    "louvain": louvain,
}


@timeit
def process_pipeline(db: Session, dataset_id: int, acquisition_ids: Sequence[int], steps: Sequence[Any]):
    """
    Process pipeline steps
    """

    dataset = dataset_service.get(db, id=dataset_id)
    cell_input = dataset.meta.get("cell")

    if not cell_input or len(acquisition_ids) == 0:
        raise HTTPException(status_code=400, detail="The dataset does not have a proper input.")

    result_create_params = ResultCreateDto(
        dataset_id=dataset_id,
        status="ready",
        name=str(datetime.utcnow()),
        pipeline=steps,
        input=acquisition_ids,
    )
    result = result_service.create(db, params=result_create_params)

    # Set ScanPy default settings
    sc.settings.autosave = False
    sc.settings.autoshow = False
    sc.settings.figdir = result.location
    sc.settings.datasetdir = result.location
    sc.settings.cachedir = result.location
    sc.settings.writedir = result.location

    adata = sc.read_h5ad(cell_input.get("location"))
    output = dict()

    # Subset observations for selected acquisitions
    adata = acquisitions_filter.process(adata, acquisition_ids=acquisition_ids, output=output)

    for step in steps:
        step_type = step.get("type")
        processor = processors.get(step_type)
        if processor:
            adata = processor.process(adata, step=step, output=output)

    # location = os.path.join(result.location, f"output{ANNDATA_FILE_EXTENSION}")
    sc.write("output", adata=adata)

    result_update_params = ResultUpdateDto(
        output=output
    )
    result = result_service.update(db, item=result, params=result_update_params)

    redis_manager.publish(
        UPDATES_CHANNEL_NAME, Message(dataset.project_id, "result_ready", jsonable_encoder(result)),
    )
