import os
import pickle
from datetime import datetime
from typing import List, Optional

import pandas as pd
from fastapi import HTTPException
from sklearn.manifold import TSNE
from sqlalchemy.orm import Session

from app.core.utils import timeit
from app.modules.dataset import crud as dataset_crud


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
    result = tsne.fit_transform(df[features].values * 2 ** 16)

    dataset = dataset_crud.get(db, id=dataset_id)
    output = dataset.output if dataset.output else {}
    tsne_output = output.get('tsne') if 'tsne' in output else []
    os.makedirs(os.path.join(dataset.location, 'tsne'), exist_ok=True)

    timestamp = str(datetime.utcnow())
    location = os.path.join(dataset.location, 'tsne', f'{timestamp}.pickle')
    with open(location, 'wb') as f:
        pickle.dump(result, f)

    tsne_output.append({
        "name": timestamp,
        "location": location,
    })

    output['tsne'] = tsne_output
    dataset_crud.update_output(db, item=dataset, output=output)
