from typing import Any, Dict

from anndata import AnnData
import numpy as np
import scanpy as sc


def process(adata: AnnData, step: Dict[str, Any]):
    """
    Transform data
    """
    mode = step.get("mode")
    if mode == "arcsinh":
        cofactor = step.get("cofactor")
        adata.X = np.arcsinh(adata.X / cofactor, out=adata.X)
    elif mode == "log1p":
        adata = sc.pp.log1p(adata, copy=True)
    return adata
