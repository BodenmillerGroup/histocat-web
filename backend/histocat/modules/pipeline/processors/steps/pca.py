from typing import Any, Dict

from anndata import AnnData
import scanpy as sc


def process(adata: AnnData, step: Dict[str, Any]):
    """
    Calculate Principal Component Analysis data
    """
    solver = step.get("svdSolver")
    adata = sc.tl.pca(adata, svd_solver=solver, copy=True)
    return adata
