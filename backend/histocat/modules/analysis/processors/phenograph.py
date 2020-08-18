import os
import pickle
from datetime import datetime
from typing import List

import pandas as pd
import phenograph
from fastapi import HTTPException
from sklearn import preprocessing
from sqlalchemy.orm import Session

from histocat.core.notifier import Message
from histocat.core.redis_manager import UPDATES_CHANNEL_NAME, redis_manager
from histocat.core.utils import timeit
from histocat.modules.dataset import service as dataset_crud


@timeit
def process_phenograph(
    db: Session,
    dataset_id: int,
    acquisition_ids: List[int],
    markers: List[str],
    clustering_algo: str,
    nearest_neighbors: int,
    jaccard: bool,
    primary_metric: str,
    min_cluster_size: int,
):
    """
    Calculate PhenoGraph data
    """

    dataset = dataset_crud.get(db, id=dataset_id)
    cell_input = dataset.input.get("cell")

    if not cell_input or len(acquisition_ids) == 0:
        raise HTTPException(status_code=400, detail="The dataset does not have a proper input.")

    df = pd.read_feather(cell_input.get("location"))
    df = df[df["AcquisitionId"].isin(acquisition_ids)]

    # Get a numpy array instead of DataFrame
    feature_values = df[markers].values

    # Normalize data
    min_max_scaler = preprocessing.MinMaxScaler()
    feature_values_scaled = min_max_scaler.fit_transform(feature_values)

    # PhenoGraph implementation
    communities, graph, Q = phenograph.cluster(
        feature_values_scaled,
        clustering_algo=clustering_algo,
        n_jobs=1,
        k=nearest_neighbors,
        jaccard=jaccard,
        primary_metric=primary_metric,
        min_cluster_size=min_cluster_size,
    )

    result = df[markers]
    result = result.assign(community=communities)
    result = result.groupby("community", as_index=False).mean()

    timestamp = str(datetime.utcnow())

    os.makedirs(os.path.join(dataset.location, "phenograph"), exist_ok=True)
    location = os.path.join(dataset.location, "phenograph", f"{timestamp}.pickle")
    with open(location, "wb") as f:
        pickle.dump(result, f)

    result = {
        "name": timestamp,
        "params": {
            "dataset_id": dataset_id,
            "acquisition_ids": acquisition_ids,
            "markers": markers,
            "clustering_algo": clustering_algo,
            "nearest_neighbors": nearest_neighbors,
            "jaccard": jaccard,
            "primary_metric": primary_metric,
            "min_cluster_size": min_cluster_size,
        },
        "location": location,
    }
    dataset_crud.update_output(db, dataset_id=dataset_id, result_type="phenograph", result=result)
    redis_manager.publish(
        UPDATES_CHANNEL_NAME, Message(dataset.experiment_id, "phenograph_result_ready", result),
    )


def get_phenograph_result(db: Session, dataset_id: int, name: str):
    """
    Read PhenoGraph result data
    """

    dataset = dataset_crud.get(db, id=dataset_id)
    phenograph_output = dataset.output.get("phenograph")

    if not phenograph_output or name not in phenograph_output:
        raise HTTPException(
            status_code=400, detail="The dataset does not have a proper PhenoGraph output.",
        )

    phenograph_result = phenograph_output.get(name)
    with open(phenograph_result.get("location"), "rb") as f:
        result = pickle.load(f)

    params = phenograph_result.get("params")
    acquisition_ids = params.get("acquisition_ids")

    output = {}
    for (columnName, columnData) in result.iteritems():
        output[columnName] = columnData.tolist()

    return output
