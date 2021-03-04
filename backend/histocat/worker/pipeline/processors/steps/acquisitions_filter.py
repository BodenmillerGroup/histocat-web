from typing import Any, Dict, Sequence

from anndata import AnnData


def process(adata: AnnData, acquisition_ids: Sequence[int], output: Dict[str, Any]):
    """
    Subset observations for selected acquisitions
    """
    result = adata[adata.obs["AcquisitionId"].isin(acquisition_ids)]
    return result
