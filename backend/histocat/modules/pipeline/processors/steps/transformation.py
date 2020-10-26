from typing import Any, Dict

import numpy as np
import scanpy as sc
from anndata import AnnData


def process(adata: AnnData, step: Dict[str, Any], output: Dict[str, Any]):
    """
    Transform data
    """
    output["transformation"] = True

    mode = step.get("mode")
    if mode == "arcsinh":
        cofactor = step.get("cofactor")
        adata.X = np.arcsinh(adata.X / cofactor, out=adata.X)
    elif mode == "log1p":
        adata = sc.pp.log1p(adata, copy=True)
    return adata
