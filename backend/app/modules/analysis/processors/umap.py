import os
from datetime import datetime
from typing import List, Optional

import numpy as np
import pandas as pd
from fastapi import HTTPException
from sklearn import preprocessing
from sqlalchemy.orm import Session
from umap import UMAP

from app.core.notifier import Message
from app.core.redis_manager import UPDATES_CHANNEL_NAME, redis_manager
from app.core.utils import timeit
from app.modules.dataset import service as dataset_crud


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

    dataset = dataset_crud.get(db, id=dataset_id)
    cell_input = dataset.input.get("cell")
    channel_map = dataset.input.get("channel_map")
    image_map = dataset.input.get("image_map")

    image_numbers = []
    for acquisition_id in acquisition_ids:
        image_number = image_map.get(str(acquisition_id))
        image_numbers.append(image_number)

    if not cell_input or not channel_map or len(image_numbers) == 0:
        raise HTTPException(status_code=400, detail="The dataset does not have a proper input.")

    df = pd.read_feather(cell_input.get("location"))
    df = df[df["ImageNumber"].isin(image_numbers)]

    features = []
    for marker in markers:
        features.append(f"Intensity_MeanIntensity_FullStack_c{channel_map[marker]}")

    # Get a numpy array instead of DataFrame
    feature_values = df[features].values

    # Normalize data
    min_max_scaler = preprocessing.MinMaxScaler()
    feature_values_scaled = min_max_scaler.fit_transform(feature_values)

    # umap-learn implementation
    umap = UMAP(
        n_components=n_components,
        n_neighbors=n_neighbors,
        min_dist=min_dist,
        metric=metric,
        verbose=6,
        random_state=42,
    )
    umap_result = umap.fit_transform(feature_values_scaled)

    timestamp = str(datetime.utcnow())

    os.makedirs(os.path.join(dataset.location, "umap"), exist_ok=True)
    location = os.path.join(dataset.location, "umap", f"{timestamp}.npy")
    np.save(location, umap_result)

    result = {
        "name": timestamp,
        "params": {
            "dataset_id": dataset_id,
            "acquisition_ids": acquisition_ids,
            "n_components": n_components,
            "n_neighbors": n_neighbors,
            "metric": metric,
            "min_dist": min_dist,
            "markers": markers,
        },
        "location": location,
    }
    dataset_crud.update_output(db, dataset_id=dataset_id, result_type="umap", result=result)
    redis_manager.publish(
        UPDATES_CHANNEL_NAME, Message(dataset.experiment_id, "umap_result_ready", result),
    )


def get_umap_result(
    db: Session, dataset_id: int, name: str, heatmap_type: Optional[str], heatmap: Optional[str],
):
    """
    Read t-SNE result data
    """

    dataset = dataset_crud.get(db, id=dataset_id)
    umap_output = dataset.output.get("umap")

    if not umap_output or name not in umap_output:
        raise HTTPException(status_code=400, detail="The dataset does not have a proper UMAP output.")

    umap_result = umap_output.get(name)
    result = np.load(umap_result.get("location"), allow_pickle=True)

    output = {
        "x": {"label": "C1", "data": result[:, 0].tolist()},
        "y": {"label": "C2", "data": result[:, 1].tolist()},
    }

    n_component = umap_result.get("params").get("n_components")
    if n_component == 3:
        output["z"] = {"label": "C3", "data": result[:, 2].tolist()}

    params = umap_result.get("params")
    acquisition_ids = params.get("acquisition_ids")
    image_map = dataset.input.get("image_map")
    cell_input = dataset.input.get("cell")

    image_numbers = []
    for acquisition_id in acquisition_ids:
        image_number = image_map.get(str(acquisition_id))
        image_numbers.append(image_number)

    df = pd.read_feather(cell_input.get("location"))
    df = df[df["ImageNumber"].isin(image_numbers)]

    if heatmap_type and heatmap:
        if heatmap_type == "channel":
            channel_map = dataset.input.get("channel_map")
            heatmap_data = df[f"Intensity_MeanIntensity_FullStack_c{channel_map[heatmap]}"] * 2 ** 16
        else:
            heatmap_data = df[heatmap]

        output["heatmap"] = {"label": heatmap, "data": heatmap_data.tolist()}
    elif len(acquisition_ids) > 1:
        image_map_inv = {v: k for k, v in image_map.items()}
        output["heatmap"] = {
            "label": "Acquisition",
            "data": [image_map_inv.get(item) for item in df["ImageNumber"]],
        }

    return output