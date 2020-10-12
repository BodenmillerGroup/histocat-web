from typing import Sequence

from anndata import AnnData


def filter_acquisitions(adata: AnnData, acquisition_ids: Sequence[int]):
    """
    Calculate Principal Component Analysis data
    """
    output = adata[adata.obs["AcquisitionId"].isin(acquisition_ids)]
    return output
