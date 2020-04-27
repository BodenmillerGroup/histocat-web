from typing import List, Optional

import pandas as pd
from fastapi import HTTPException
from sklearn import preprocessing
from sklearn.decomposition import PCA
from sqlalchemy.orm import Session

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

    pca = PCA(n_components=n_components)

    # Get a numpy array instead of DataFrame
    feature_values = df[features].values

    # Normalize data
    min_max_scaler = preprocessing.MinMaxScaler()
    feature_values_scaled = min_max_scaler.fit_transform(feature_values)

    # Run PCA
    result = pca.fit_transform(feature_values_scaled)

    cell_ids = df["ImageNumber"].astype(str) + "_" + df["ObjectNumber"].astype(str)
    output = {
        "cell_ids": cell_ids.tolist(),
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
    elif len(acquisition_ids) > 1:
        image_map_inv = {v: k for k, v in image_map.items()}
        output["heatmap"] = {
            "label": "Acquisition",
            "data": [image_map_inv.get(item) for item in df["ImageNumber"].tolist()],
        }

    return output
