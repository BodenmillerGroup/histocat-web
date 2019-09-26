import os
from datetime import datetime
import pickle
from typing import List, Optional

import numpy as np
import pandas as pd
from scipy import sparse
from sklearn import preprocessing
from fastapi import HTTPException
from sqlalchemy.orm import Session
import phenograph

from app.core.notifier import Message
from app.core.redis_manager import redis_manager, UPDATES_CHANNEL_NAME
from app.core.utils import timeit
from app.modules.dataset import crud as dataset_crud


@timeit
def process_phenograph(
    db: Session,
    dataset_id: int,
    acquisition_ids: List[int],
    markers: List[str],
):
    """
    Calculate PhenoGraph data
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
        raise HTTPException(
            status_code=400,
            detail="The dataset does not have a proper input.",
        )

    df = pd.read_feather(cell_input.get("location"))
    df = df[df["ImageNumber"].isin(image_numbers)]

    features = []
    for marker in markers:
        features.append(f'Intensity_MeanIntensity_FullStack_c{channel_map[marker]}')

    # Get a numpy array instead of DataFrame
    feature_values = df[features].values

    # Normalize data
    min_max_scaler = preprocessing.MinMaxScaler()
    feature_values_scaled = min_max_scaler.fit_transform(feature_values)

    # PhenoGraph implementation
    communities, graph, Q = phenograph.cluster(feature_values_scaled, n_jobs=1)

    timestamp = str(datetime.utcnow())

    os.makedirs(os.path.join(dataset.location, 'phenograph'), exist_ok=True)

    communities_location = os.path.join(dataset.location, 'phenograph', f'{timestamp}_communities.npy')
    np.save(communities_location, communities)

    graph_location = os.path.join(dataset.location, 'phenograph', f'{timestamp}_graph.npz')
    sparse.save_npz(graph_location, graph)

    q_location = os.path.join(dataset.location, 'phenograph', f'{timestamp}_Q.pickle')
    with open(q_location, "wb") as f:
        pickle.dump(Q, f)

    result = {
        "name": timestamp,
        "params": {
            "dataset_id": dataset_id,
            "acquisition_ids": acquisition_ids,
            "markers": markers,
        },
        "communities_location": communities_location,
        "graph_location": graph_location,
        "q_location": q_location,
    }
    dataset_crud.update_output(db, dataset_id=dataset_id, result_type='phenograph', result=result)
    redis_manager.publish(UPDATES_CHANNEL_NAME, Message(dataset.experiment_id, "phenograph_result_ready", result))


def get_phenograph_result(
    db: Session,
    dataset_id: int,
    name: str,
    heatmap_type: Optional[str],
    heatmap: Optional[str],
):
    """
    Read PhenoGraph result data
    """

    dataset = dataset_crud.get(db, id=dataset_id)
    phenograph_output = dataset.output.get("phenograph")

    if not phenograph_output or name not in phenograph_output:
        raise HTTPException(
            status_code=400,
            detail="The dataset does not have a proper PhenoGraph output.",
        )

    phenograph_result = phenograph_output.get(name)
    communities = np.load(phenograph_result.get("communities_location"), allow_pickle=True)
    graph = sparse.load_npz(phenograph_result.get("graph_location"))
    with open(phenograph_result.get("q_location"), "rb") as f:
        Q = pickle.load(f)

    output = {
        "communities": {
            "label": "communities",
            "data": communities
        },
        "graph": {
            "label": "graph",
            "data": graph
        },
        "Q": {
            "label": "Q",
            "data": Q
        }
    }

    # params = umap_result.get("params")
    # acquisition_ids = params.get("acquisition_ids")
    # image_map = dataset.input.get("image_map")
    # cell_input = dataset.input.get("cell")
    #
    # image_numbers = []
    # for acquisition_id in acquisition_ids:
    #     image_number = image_map.get(str(acquisition_id))
    #     image_numbers.append(image_number)
    #
    # df = pd.read_feather(cell_input.get("location"))
    # df = df[df["ImageNumber"].isin(image_numbers)]
    #
    # if heatmap_type and heatmap:
    #     if heatmap_type == "channel":
    #         channel_map = dataset.input.get("channel_map")
    #         heatmap_data = df[f'Intensity_MeanIntensity_FullStack_c{channel_map[heatmap]}'] * 2 ** 16
    #     else:
    #         heatmap_data = df[heatmap]
    #
    #     output["heatmap"] = {
    #         "label": heatmap,
    #         "data": heatmap_data
    #     }
    # elif len(acquisition_ids) > 1:
    #     image_map_inv = {v: k for k, v in image_map.items()}
    #     output["heatmap"] = {
    #         "label": "Acquisition",
    #         "data": [image_map_inv.get(item) for item in df["ImageNumber"]]
    #     }

    return output
