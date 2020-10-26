from typing import Any, Dict

import scanpy as sc
from anndata import AnnData


def process(adata: AnnData, step: Dict[str, Any], output: Dict[str, Any]):
    """
    Compute a neighborhood graph of observations
    """
    output["neighbors"] = True

    n_neighbors = step.get("nNeighbors")
    knn = step.get("knn")
    method = step.get("method")
    metric = step.get("metric")
    result = sc.pp.neighbors(adata, n_neighbors=n_neighbors, knn=knn, method=method, metric=metric, copy=True)
    return result
