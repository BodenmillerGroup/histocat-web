from typing import List, Optional

import pandas as pd
from fastapi import HTTPException
from sklearn.decomposition import PCA
from sqlalchemy.orm import Session

from app.modules.dataset import crud as dataset_crud


def process_pca(
    db: Session,
    dataset_id: int,
    acquisition_id: int,
    n_components: int,
    markers: List[str],
    heatmap: Optional[str],
):
    """
    Calculate Principal Component Analysis data
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

    pca = PCA(n_components=n_components)
    result = pca.fit_transform(df[features].values * 2 ** 16)

    output = {
        "x": {
            "label": "PC1",
            "data": result[:, 0]
        },
        "y": {
            "label": "PC2",
            "data": result[:, 1],
        },
    }

    if n_components == 3:
        output["z"] = {
            "label": "PC3",
            "data": result[:, 2],
        }

    if heatmap:
        output["heatmap"] = {
            "label": heatmap,
            "data": df[heatmap]
        }

    return output