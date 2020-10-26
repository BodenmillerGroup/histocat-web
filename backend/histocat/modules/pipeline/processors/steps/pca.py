from typing import Any, Dict

import scanpy as sc
from anndata import AnnData


def process(adata: AnnData, step: Dict[str, Any], output: Dict[str, Any]):
    """
    Calculate Principal Component Analysis data
    """
    output["pca"] = True

    solver = step.get("svdSolver")
    result = sc.tl.pca(adata, n_comps=2, svd_solver=solver, copy=True)
    return result
