import os
import pickle
from datetime import datetime
from typing import List, Optional

import pandas as pd
import numpy as np
from fastapi import HTTPException
from sklearn.manifold import TSNE
from sqlalchemy.orm import Session

from app.core.notifier import Message
from app.core.utils import timeit
from app.modules.dataset import crud as dataset_crud
from app.core.redis_manager import redis_manager, UPDATES_CHANNEL_NAME


@timeit
def process_tsne(
    db: Session,
    dataset_id: int,
    acquisition_id: int,
    n_components: int,
    markers: List[str],
    heatmap: Optional[str],
):
    """
    Calculate t-Distributed Stochastic Neighbor Embedding data
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

    tsne = TSNE(n_components=n_components)
    tsne_result = tsne.fit_transform(df[features].values * 2 ** 16)

    timestamp = str(datetime.utcnow())

    os.makedirs(os.path.join(dataset.location, 'tsne'), exist_ok=True)
    location = os.path.join(dataset.location, 'tsne', f'{timestamp}.npy')
    np.save(location, tsne_result)

    result = {
        "name": timestamp,
        "params": {
            "dataset_id": dataset_id,
            "acquisition_id": acquisition_id,
            "n_components": n_components,
            "markers": markers,
            "heatmap": heatmap,
        },
        "location": location,
    }
    dataset_crud.update_output(db, dataset_id=dataset_id, result=result)
    redis_manager.publish(UPDATES_CHANNEL_NAME, Message(dataset.experiment_id, "tsne_result_ready", result))


def get_tsne_result(
    db: Session,
    dataset_id: int,
    name: str,
):
    """
    Read t-SNE result data
    """

    dataset = dataset_crud.get(db, id=dataset_id)
    tsne_output = dataset.output.get("tsne")

    if not tsne_output or name not in tsne_output:
        raise HTTPException(
            status_code=400,
            detail="The dataset does not have a proper t-SNE output.",
        )

    tsne_result = tsne_output.get(name)
    result = np.load(tsne_result.get("location"), allow_pickle=True)

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

    n_component = tsne_result.get("params").get("n_components")
    if n_component == 3:
        output["z"] = {
            "label": "C3",
            "data": result[:, 2].tolist()
        }

    params = tsne_result.get("params")
    heatmap = params.get("heatmap")
    if heatmap:
        cell_input = dataset.input.get("cell")
        image_map = dataset.input.get("image_map")
        acquisition_id = params.get("acquisition_id")
        image_number = image_map.get(str(acquisition_id))
        df = pd.read_feather(cell_input.get("location"))
        df = df[df["ImageNumber"] == image_number]
        output["heatmap"] = {
            "label": heatmap,
            "data": df[heatmap]
        }

    return output
