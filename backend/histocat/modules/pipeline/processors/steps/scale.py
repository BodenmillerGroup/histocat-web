from typing import Any, Dict

import scanpy as sc
from anndata import AnnData


def process(adata: AnnData, step: Dict[str, Any], output: Dict[str, Any]):
    """
    Scale data
    """
    output["scale"] = True

    zero_center = step.get("zeroCenter")
    max_value = step.get("maxValue")
    result = sc.pp.scale(adata, zero_center=zero_center, max_value=max_value, copy=True)
    return result
