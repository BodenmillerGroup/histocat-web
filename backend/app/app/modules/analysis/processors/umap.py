import os
from datetime import datetime
from typing import List, Optional

import numpy as np
import pandas as pd
from fastapi import HTTPException
from sqlalchemy.orm import Session
from umap import UMAP

from app.core.notifier import Message
from app.core.redis_manager import redis_manager, UPDATES_CHANNEL_NAME
from app.core.utils import timeit
from app.modules.dataset import crud as dataset_crud


@timeit
def process_umap(
    db: Session,
    dataset_id: int,
    acquisition_id: int,
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
    image_number = image_map.get(str(acquisition_id))
    if not cell_input or not image_number or not channel_map:
        raise HTTPException(
            status_code=400,
            detail="The dataset does not have a proper input.",
        )

    df = pd.read_feather(cell_input.get("location"))
    df = df[df["ImageNumber"] == image_number]

    features = []
    for marker in markers:
        features.append(f'Intensity_MeanIntensity_FullStack_c{channel_map[marker]}')

    # umap-learn implementation
    umap = UMAP(n_components=n_components, n_neighbors=n_neighbors, min_dist=min_dist, metric=metric, verbose=6, random_state=42)
    umap_result = umap.fit_transform(df[features].values * 2 ** 16)

    timestamp = str(datetime.utcnow())

    os.makedirs(os.path.join(dataset.location, 'umap'), exist_ok=True)
    location = os.path.join(dataset.location, 'umap', f'{timestamp}.npy')
    np.save(location, umap_result)

    result = {
        "name": timestamp,
        "params": {
            "dataset_id": dataset_id,
            "acquisition_id": acquisition_id,
            "n_components": n_components,
            "n_neighbors": n_neighbors,
            "metric": metric,
            "min_dist": min_dist,
            "markers": markers,
        },
        "location": location,
    }
    dataset_crud.update_output(db, dataset_id=dataset_id, result_type='umap', result=result)
    redis_manager.publish(UPDATES_CHANNEL_NAME, Message(dataset.experiment_id, "umap_result_ready", result))


def get_umap_result(
    db: Session,
    dataset_id: int,
    name: str,
    heatmap_type: Optional[str],
    heatmap: Optional[str],
):
    """
    Read t-SNE result data
    """

    dataset = dataset_crud.get(db, id=dataset_id)
    umap_output = dataset.output.get("umap")

    if not umap_output or name not in umap_output:
        raise HTTPException(
            status_code=400,
            detail="The dataset does not have a proper UMAP output.",
        )

    umap_result = umap_output.get(name)
    result = np.load(umap_result.get("location"), allow_pickle=True)

    output = {
        "x": {
            "label": "C1",
            "data": result[:, 0].tolist()
        },
        "y": {
            "label": "C2",
            "data": result[:, 1].tolist()
        },
    }

    n_component = umap_result.get("params").get("n_components")
    if n_component == 3:
        output["z"] = {
            "label": "C3",
            "data": result[:, 2].tolist()
        }

    if heatmap_type and heatmap:
        params = umap_result.get("params")
        cell_input = dataset.input.get("cell")
        image_map = dataset.input.get("image_map")
        acquisition_id = params.get("acquisition_id")
        image_number = image_map.get(str(acquisition_id))

        df = pd.read_feather(cell_input.get("location"))
        df = df[df["ImageNumber"] == image_number]

        if heatmap_type == "channel":
            channel_map = dataset.input.get("channel_map")
            heatmap_data = df[f'Intensity_MeanIntensity_FullStack_c{channel_map[heatmap]}'] * 2 ** 16
        else:
            heatmap_data = df[heatmap]

        output["heatmap"] = {
            "label": heatmap,
            "data": heatmap_data
        }

    return output
