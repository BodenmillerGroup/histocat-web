from typing import List, Optional

import pandas as pd
from fastapi import HTTPException
from sklearn.decomposition import PCA
from sqlalchemy.orm import Session

from app.core.image import get_heatmap_colors
from app.modules.dataset import crud as dataset_crud


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
        raise HTTPException(
            status_code=400,
            detail="The dataset does not have a proper input.",
        )

    df = pd.read_feather(cell_input.get("location"))
    df = df[df["ImageNumber"].isin(image_numbers)]

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

    if heatmap_type and heatmap:
        if heatmap_type == "channel":
            channel_map = dataset.input.get("channel_map")
            heatmap_data = df[f'Intensity_MeanIntensity_FullStack_c{channel_map[heatmap]}'] * 2 ** 16
        else:
            heatmap_data = df[heatmap]

        output["heatmap"] = {
            "label": heatmap,
            "data": heatmap_data
        }
    elif len(acquisition_ids) > 1:
        output["heatmap"] = {
            "label": "Acquisition",
            "data": get_heatmap_colors(df["ImageNumber"], True)
        }

    return output
