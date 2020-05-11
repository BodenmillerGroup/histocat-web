import os
import pickle
from datetime import datetime
from typing import List, Optional

import numpy as np
import pandas as pd
from fastapi import HTTPException
from sklearn import preprocessing
from sklearn.manifold import TSNE

from sqlalchemy.orm import Session

from histocat.core.image import normalize_embedding
from histocat.core.notifier import Message
from histocat.core.redis_manager import UPDATES_CHANNEL_NAME, redis_manager
from histocat.core.utils import timeit
from histocat.modules.dataset import service as dataset_crud


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

    dataset = dataset_crud.get(db, id=dataset_id)
    cell_input = dataset.input.get("cell")
    channel_map = dataset.input.get("channel_map")

    if not cell_input or not channel_map or len(acquisition_ids) == 0:
        raise HTTPException(status_code=400, detail="The dataset does not have a proper input.")

    df = pd.read_feather(cell_input.get("location"))
    df = df[df["acquisition_id"].isin(acquisition_ids)]

    features = []
    for marker in markers:
        features.append(f"Intensity_MeanIntensity_FullStack_c{channel_map[marker]}")

    # Get a numpy array instead of DataFrame
    feature_values = df[features].values

    # Normalize data
    feature_values = np.arcsinh(feature_values / 5, out=feature_values)

    min_max_scaler = preprocessing.MinMaxScaler()
    feature_values_scaled = min_max_scaler.fit_transform(feature_values)

    # scikit-learn implementation
    tsne = TSNE(
        n_components=n_components,
        perplexity=perplexity,
        learning_rate=learning_rate,
        n_iter=iterations,
        verbose=6,
        random_state=42,
        init=init,
    )
    tsne_result = tsne.fit_transform(feature_values_scaled)
    cell_ids = df["acquisition_id"].astype(str) + "_" + df["ObjectNumber"].astype(str)

    timestamp = str(datetime.utcnow())

    os.makedirs(os.path.join(dataset.location, "tsne"), exist_ok=True)
    location = os.path.join(dataset.location, "tsne", f"{timestamp}.pickle")

    with open(location, "wb") as f:
        pickle.dump({"cell_ids": cell_ids, "tsne_result": tsne_result}, f, pickle.HIGHEST_PROTOCOL)

    result = {
        "name": timestamp,
        "params": {
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
        "location": location,
    }
    dataset_crud.update_output(db, dataset_id=dataset_id, result_type="tsne", result=result)
    redis_manager.publish(
        UPDATES_CHANNEL_NAME, Message(dataset.experiment_id, "tsne_result_ready", result),
    )


def get_tsne_result(
    db: Session, dataset_id: int, name: str, heatmap_type: Optional[str], heatmap: Optional[str],
):
    """
    Read t-SNE result data
    """

    dataset = dataset_crud.get(db, id=dataset_id)
    tsne_output = dataset.output.get("tsne")

    if not tsne_output or name not in tsne_output:
        raise HTTPException(status_code=400, detail="The dataset does not have a proper t-SNE output.")

    tsne_result = tsne_output.get(name)

    with open(tsne_result.get("location"), "rb") as f:
        r = pickle.load(f)

    cell_ids = r.get("cell_ids")
    result = r.get("tsne_result")

    result = normalize_embedding(result)

    output = {
        "cell_ids": cell_ids.tolist(),
        "x": {"label": "C1", "data": result[:, 0].tolist()},
        "y": {"label": "C2", "data": result[:, 1].tolist()},
    }

    n_component = tsne_result.get("params").get("n_components")
    if n_component == 3:
        output["z"] = {"label": "C3", "data": result[:, 2].tolist()}

    params = tsne_result.get("params")
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
            heatmap_data = (df[f"Intensity_MeanIntensity_FullStack_c{channel_map[heatmap]}"] * 2 ** 16).tolist()
        else:
            heatmap_data = df[heatmap].tolist()

        output["heatmap"] = {"label": heatmap, "data": heatmap_data}
    elif len(acquisition_ids) > 1:
        image_map_inv = {v: k for k, v in image_map.items()}
        output["heatmap"] = {
            "label": "Acquisition",
            "data": [image_map_inv.get(item) for item in df["ImageNumber"]],
        }

    return output
