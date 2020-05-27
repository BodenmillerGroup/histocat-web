from typing import List, Optional

import pandas as pd
from fastapi import HTTPException
from sklearn import preprocessing
from sklearn.decomposition import PCA
from sqlalchemy.orm import Session

from histocat.core.image import normalize_embedding
from histocat.modules.dataset import service as dataset_crud


def process_pca(
    db: Session,
    dataset_id: int,
    acquisition_ids: List[int],
    n_components: int,
    markers: List[str],
    heatmap_type: Optional[str],
    heatmap: Optional[str],
):
    """
    Calculate Principal Component Analysis data
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

    pca = PCA(n_components=n_components)

    # Get a numpy array instead of DataFrame
    feature_values = df[features].values

    # Normalize data
    min_max_scaler = preprocessing.MinMaxScaler()
    feature_values_scaled = min_max_scaler.fit_transform(feature_values)

    # Run PCA
    result = pca.fit_transform(feature_values_scaled)
    result = normalize_embedding(result)

    output = {
        "acquisitionIds": df["acquisition_id"].tolist(),
        "cellIds": df["ObjectNumber"].tolist(),
        "x": {"label": "PC1", "data": result[:, 0].tolist()},
        "y": {"label": "PC2", "data": result[:, 1].tolist()},
    }

    if n_components == 3:
        output["z"] = {"label": "PC3", "data": result[:, 2].tolist()}

    if heatmap_type and heatmap:
        if heatmap_type == "channel":
            channel_map = dataset.input.get("channel_map")
            heatmap_data = df[f"Intensity_MeanIntensity_FullStack_c{channel_map[heatmap]}"] * 2 ** 16
        else:
            heatmap_data = df[heatmap]

        output["heatmap"] = {"label": heatmap, "data": heatmap_data.tolist()}

    return output
