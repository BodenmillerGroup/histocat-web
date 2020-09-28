import os
import pickle
from datetime import datetime
from typing import List, Optional

import numpy as np
import pandas as pd
import scanpy as sc
from fastapi import HTTPException
from fastapi.encoders import jsonable_encoder
from sklearn import preprocessing
from sklearn.manifold import TSNE
from sqlalchemy.orm import Session

from histocat.core.image import normalize_embedding
from histocat.core.notifier import Message
from histocat.core.redis_manager import UPDATES_CHANNEL_NAME, redis_manager
from histocat.core.utils import timeit
from histocat.modules.dataset import service as dataset_service
from histocat.modules.result import service as result_service
from histocat.modules.result.dto import ResultCreateDto


@timeit
def process_tsne(
    db: Session,
    dataset_id: int,
    acquisition_ids: List[int],
    n_components: int,
    perplexity: int,
    learning_rate: int,
    iterations: int,
    theta: float,
    init: str,
    markers: List[str],
):
    """
    Calculate t-Distributed Stochastic Neighbor Embedding data
    """

    dataset = dataset_service.get(db, id=dataset_id)
    cell_input = dataset.meta.get("cell")

    if not cell_input or len(acquisition_ids) == 0:
        raise HTTPException(status_code=400, detail="The dataset does not have a proper input.")

    adata = sc.read_h5ad(cell_input.get("location"))

    # Subset observations for selected acquisitions
    adata = adata[adata.obs["AcquisitionId"].isin(acquisition_ids)]

    # Subset selected channels
    feature_values = adata[:, adata.var.index.isin(markers)].layers["expr"]

    # scikit-learn implementation
    tsne = TSNE(
        n_components=n_components,
        perplexity=perplexity,
        learning_rate=learning_rate,
        n_iter=iterations,
        verbose=6,
        random_state=77,
        init=init,
    )
    tsne_result = tsne.fit_transform(feature_values)
    acquisitionIds = adata.obs["AcquisitionId"].tolist()
    cellIds = adata.obs["CellId"].tolist()
    objectNumbers = adata.obs["ObjectNumber"].tolist()

    result_create_params = ResultCreateDto(
        dataset_id=dataset_id,
        type="tsne",
        status="ready",
        name=str(datetime.utcnow()),
        params={
            "dataset_id": dataset_id,
            "acquisition_ids": acquisition_ids,
            "n_components": n_components,
            "perplexity": perplexity,
            "learning_rate": learning_rate,
            "iterations": iterations,
            "theta": theta,
            "init": init,
            "markers": markers,
        },
    )
    result = result_service.create(db, params=result_create_params)

    location = os.path.join(result.location, "data.pickle")
    with open(location, "wb") as f:
        pickle.dump(
            {
                "acquisitionIds": acquisitionIds,
                "cellIds": cellIds,
                "objectNumbers": objectNumbers,
                "tsne_result": tsne_result,
            },
            f,
            pickle.HIGHEST_PROTOCOL,
        )

    redis_manager.publish(
        UPDATES_CHANNEL_NAME, Message(dataset.project_id, "tsne_result_ready", jsonable_encoder(result)),
    )


def get_tsne_result(
    db: Session, result_id: int, heatmap_type: Optional[str], heatmap: Optional[str],
):
    """
    Read t-SNE result data
    """

    result = result_service.get(db, id=result_id)
    if not result:
        raise HTTPException(status_code=404, detail="t-SNE result not found.")

    location = os.path.join(result.location, "data.pickle")
    with open(location, "rb") as f:
        r = pickle.load(f)

    acquisitionIds = r.get("acquisitionIds")
    cellIds = r.get("cellIds")
    objectNumbers = r.get("objectNumbers")
    tsne_result = r.get("tsne_result")

    output = {
        "acquisitionIds": acquisitionIds,
        "cellIds": cellIds,
        "objectNumbers": objectNumbers,
        "x": {"label": "C1", "data": tsne_result[:, 0].tolist()},
        "y": {"label": "C2", "data": tsne_result[:, 1].tolist()},
    }

    params = result.params
    n_component = params.get("n_components")
    if n_component == 3:
        output["z"] = {"label": "C3", "data": tsne_result[:, 2].tolist()}

    if heatmap_type and heatmap:
        acquisition_ids = params.get("acquisition_ids")
        cell_input = result.dataset.meta.get("cell")

        adata = sc.read_h5ad(cell_input.get("location"))

        # Subset observations for selected acquisitions
        adata = adata[adata.obs["AcquisitionId"].isin(acquisition_ids)]

        heatmap_values = adata.layers["expr"][:, adata.var.index == heatmap]
        output["heatmap"] = {"label": heatmap, "data": heatmap_values[:, 0].tolist()}

    return output
