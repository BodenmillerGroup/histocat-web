from typing import Any, Dict

import scanpy as sc
from anndata import AnnData


def process(adata: AnnData, step: Dict[str, Any], output: Dict[str, Any]):
    """
    Compute a neighborhood graph of observations
    """
    output["neighbors"] = True

    n_neighbors = step.get("nNeighbors")
    metric = step.get("metric")
    random_state = step.get("randomState")

    result = sc.pp.neighbors(adata, n_neighbors=n_neighbors, knn=True, method="umap", random_state=random_state, metric=metric, copy=True)
    return result
