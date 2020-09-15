import os
import pickle
from datetime import datetime
from typing import List

import pandas as pd
import phenograph
from fastapi import HTTPException
from fastapi.encoders import jsonable_encoder
from sklearn import preprocessing
from sqlalchemy.orm import Session

from histocat.core.notifier import Message
from histocat.core.redis_manager import UPDATES_CHANNEL_NAME, redis_manager
from histocat.core.utils import timeit
from histocat.modules.dataset import service as dataset_service
from histocat.modules.result import service as result_service
from histocat.modules.result.dto import ResultCreateDto


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

    dataset = dataset_service.get(db, id=dataset_id)
    cell_input = dataset.meta.get("cell")

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

    phenograph_result = df[markers]
    phenograph_result = phenograph_result.assign(community=communities)
    phenograph_result = phenograph_result.groupby("community", as_index=False).mean()

    result_create_params = ResultCreateDto(
        dataset_id=dataset_id,
        type="phenograph",
        status="ready",
        name=str(datetime.utcnow()),
        params={
            "dataset_id": dataset_id,
            "acquisition_ids": acquisition_ids,
            "markers": markers,
            "clustering_algo": clustering_algo,
            "nearest_neighbors": nearest_neighbors,
            "jaccard": jaccard,
            "primary_metric": primary_metric,
            "min_cluster_size": min_cluster_size,
        },
    )
    result = result_service.create(db, params=result_create_params)

    location = os.path.join(result.location, "data.pickle")
    with open(location, "wb") as f:
        pickle.dump(phenograph_result, f)

    redis_manager.publish(
        UPDATES_CHANNEL_NAME, Message(dataset.experiment_id, "phenograph_result_ready", jsonable_encoder(result)),
    )


def get_phenograph_result(db: Session, result_id: int):
    """
    Read PhenoGraph result data
    """

    result = result_service.get(db, id=result_id)
    if not result:
        raise HTTPException(status_code=404, detail="PhenoGraph result not found.")

    location = os.path.join(result.location, "data.pickle")
    with open(location, "rb") as f:
        r = pickle.load(f)

    output = {}
    for (columnName, columnData) in r.iteritems():
        output[columnName] = columnData.tolist()

    return output
