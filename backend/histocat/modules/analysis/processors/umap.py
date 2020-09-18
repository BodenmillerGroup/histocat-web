import os
import pickle
from datetime import datetime
from typing import List, Optional

import numpy as np
import pandas as pd
from fastapi import HTTPException
from fastapi.encoders import jsonable_encoder
from sklearn import preprocessing
from sqlalchemy.orm import Session
from umap import UMAP

from histocat.core.image import normalize_embedding
from histocat.core.notifier import Message
from histocat.core.redis_manager import UPDATES_CHANNEL_NAME, redis_manager
from histocat.core.utils import timeit
from histocat.modules.dataset import service as dataset_service
from histocat.modules.result import service as result_service
from histocat.modules.result.dto import ResultCreateDto


@timeit
def process_umap(
    db: Session,
    dataset_id: int,
    acquisition_ids: List[int],
    n_components: int,
    n_neighbors: int,
    metric: str,
    min_dist: float,
    markers: List[str],
):
    """
    Calculate Uniform Manifold Approximation and Projection data
    """

    dataset = dataset_service.get(db, id=dataset_id)
    cell_input = dataset.meta.get("cell")

    if not cell_input or len(acquisition_ids) == 0:
        raise HTTPException(status_code=400, detail="The dataset does not have a proper input.")

    df = pd.read_feather(cell_input.get("location"))
    df = df[df["AcquisitionId"].isin(acquisition_ids)]

    # Get a numpy array instead of DataFrame
    feature_values = df[markers].values

    # Normalize data
    feature_values = np.arcsinh(feature_values / 5, out=feature_values)

    min_max_scaler = preprocessing.MinMaxScaler()
    feature_values_scaled = min_max_scaler.fit_transform(feature_values)

    # umap-learn implementation
    umap = UMAP(
        n_components=n_components,
        n_neighbors=n_neighbors,
        min_dist=min_dist,
        metric=metric,
        verbose=6,
        random_state=77,
    )
    umap_result = umap.fit_transform(feature_values_scaled)
    acquisitionIds = df["AcquisitionId"]
    cellIds = df["CellId"]
    objectNumbers = df["ObjectNumber"]

    result_create_params = ResultCreateDto(
        dataset_id=dataset_id,
        type="umap",
        status="ready",
        name=str(datetime.utcnow()),
        params={
            "dataset_id": dataset_id,
            "acquisition_ids": acquisition_ids,
            "n_components": n_components,
            "n_neighbors": n_neighbors,
            "metric": metric,
            "min_dist": min_dist,
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
                "umap_result": umap_result,
            },
            f,
            pickle.HIGHEST_PROTOCOL,
        )

    redis_manager.publish(
        UPDATES_CHANNEL_NAME, Message(dataset.project_id, "umap_result_ready", jsonable_encoder(result)),
    )


def get_umap_result(
    db: Session, result_id: int, heatmap_type: Optional[str], heatmap: Optional[str],
):
    """
    Read t-SNE result data
    """

    result = result_service.get(db, id=result_id)
    if not result:
        raise HTTPException(status_code=404, detail="UMAP result not found.")

    location = os.path.join(result.location, "data.pickle")
    with open(location, "rb") as f:
        r = pickle.load(f)

    acquisitionIds = r.get("acquisitionIds")
    cellIds = r.get("cellIds")
    objectNumbers = r.get("objectNumbers")
    umap_result = r.get("umap_result")

    umap_result = normalize_embedding(umap_result)

    output = {
        "acquisitionIds": acquisitionIds.tolist(),
        "cellIds": cellIds.tolist(),
        "objectNumbers": objectNumbers.tolist(),
        "x": {"label": "C1", "data": umap_result[:, 0].tolist()},
        "y": {"label": "C2", "data": umap_result[:, 1].tolist()},
    }

    params = result.params
    n_component = params.get("n_components")
    if n_component == 3:
        output["z"] = {"label": "C3", "data": umap_result[:, 2].tolist()}

    if heatmap_type and heatmap:
        acquisition_ids = params.get("acquisition_ids")
        cell_input = result.dataset.meta.get("cell")

        df = pd.read_feather(cell_input.get("location"))
        df = df[df["AcquisitionId"].isin(acquisition_ids)]

        output["heatmap"] = {"label": heatmap, "data": df[heatmap].tolist()}

    return output
