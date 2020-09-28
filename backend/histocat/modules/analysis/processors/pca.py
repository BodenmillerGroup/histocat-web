from typing import List, Optional

import scanpy as sc
from fastapi import HTTPException
from sqlalchemy.orm import Session

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

    adata = sc.read_h5ad(cell_input.get("location"))

    # Subset observations for selected acquisitions
    adata = adata[adata.obs["AcquisitionId"].isin(acquisition_ids)]

    # Subset selected channels
    feature_values = adata[:, adata.var.index.isin(markers)].layers["expr"]

    result = sc.tl.pca(feature_values, n_comps=n_components, copy=True)

    output = {
        "acquisitionIds": adata.obs["AcquisitionId"].tolist(),
        "cellIds": adata.obs["CellId"].tolist(),
        "objectNumbers": adata.obs["ObjectNumber"].tolist(),
        "x": {"label": "PC1", "data": result[:, 0].tolist()},
        "y": {"label": "PC2", "data": result[:, 1].tolist()},
    }

    if n_components == 3:
        output["z"] = {"label": "PC3", "data": result[:, 2].tolist()}

    if heatmap_type and heatmap:
        heatmap_values = adata.layers["expr"][:, adata.var.index == heatmap]
        output["heatmap"] = {"label": heatmap, "data": heatmap_values[:, 0].tolist()}

    return output
