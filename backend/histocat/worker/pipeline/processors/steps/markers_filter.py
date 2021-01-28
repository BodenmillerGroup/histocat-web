from typing import Any, Dict

from anndata import AnnData


def process(adata: AnnData, step: Dict[str, Any], output: Dict[str, Any]):
    """
    Subset variables for selected markers
    """
    output["markersFilter"] = True

    markers = step.get("markers")
    result = adata[:, adata.var.index.isin(markers)]
    return result
