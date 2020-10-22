from typing import Any, Dict

from anndata import AnnData
import scanpy as sc


def process(adata: AnnData, step: Dict[str, Any], output: Dict[str, Any]):
    """
    Scale data
    """
    output["scale"] = True

    max_value = step.get("maxValue")
    result = sc.pp.scale(adata, max_value=max_value, copy=True)
    return result
