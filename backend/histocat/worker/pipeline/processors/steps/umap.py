from typing import Any, Dict

import scanpy as sc
from anndata import AnnData


def process(adata: AnnData, step: Dict[str, Any], output: Dict[str, Any]):
    """Calculate Uniform Manifold Approximation and Projection data"""
    output["umap"] = True

    min_dist = step.get("minDist")
    spread = step.get("spread")
    result = sc.tl.umap(adata, n_components=2, min_dist=min_dist, spread=spread, copy=True)
    return result
