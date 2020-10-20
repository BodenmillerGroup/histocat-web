from typing import Any, Dict

from anndata import AnnData


def process(adata: AnnData, step: Dict[str, Any]):
    """
    Subset variables for selected markers
    """
    markers = step.get("markers")
    output = adata[:, adata.var.index.isin(markers)]
    return output
