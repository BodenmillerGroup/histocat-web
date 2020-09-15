from typing import List, Optional

import pandas as pd
from fastapi import HTTPException
from sklearn import preprocessing
from sklearn.decomposition import PCA
from sqlalchemy.orm import Session

from histocat.core.image import normalize_embedding
from histocat.modules.dataset import service


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

    dataset = service.get(db, id=dataset_id)
    cell_input = dataset.meta.get("cell")

    if not cell_input or len(acquisition_ids) == 0:
        raise HTTPException(status_code=400, detail="The dataset does not have a proper input.")

    df = pd.read_feather(cell_input.get("location"))
    df = df[df["AcquisitionId"].isin(acquisition_ids)]

    pca = PCA(n_components=n_components)

    # Get a numpy array instead of DataFrame
    feature_values = df[markers].values

    # Normalize data
    min_max_scaler = preprocessing.MinMaxScaler()
    feature_values_scaled = min_max_scaler.fit_transform(feature_values)

    # Run PCA
    result = pca.fit_transform(feature_values_scaled)
    result = normalize_embedding(result)

    output = {
        "acquisitionIds": df["AcquisitionId"].tolist(),
        "cellIds": df["CellId"].tolist(),
        "objectNumbers": df["ObjectNumber"].tolist(),
        "x": {"label": "PC1", "data": result[:, 0].tolist()},
        "y": {"label": "PC2", "data": result[:, 1].tolist()},
    }

    if n_components == 3:
        output["z"] = {"label": "PC3", "data": result[:, 2].tolist()}

    if heatmap_type and heatmap:
        output["heatmap"] = {"label": heatmap, "data": df[heatmap].tolist()}

    return output
