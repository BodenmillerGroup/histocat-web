from typing import Sequence

from anndata import AnnData


def process(adata: AnnData, acquisition_ids: Sequence[int]):
    """
    Subset observations for selected acquisitions
    """
    output = adata[adata.obs["AcquisitionId"].isin(acquisition_ids)]
    return output
