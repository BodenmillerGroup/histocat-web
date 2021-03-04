from typing import Any, Dict

import scanpy as sc
from anndata import AnnData


def process(adata: AnnData, step: Dict[str, Any], output: Dict[str, Any]):
    """Cluster cells using the Leiden algorithm"""
    output["leiden"] = True

    resolution = step.get("resolution")
    directed = step.get("directed")
    use_weights = step.get("useWeights")
    result = sc.tl.leiden(adata, resolution=resolution, directed=directed, use_weights=use_weights, copy=True)
    sc.pl.stacked_violin(result, var_names=result.var_names, groupby=["leiden"], save="leiden.png", dendrogram=True)
    sc.pl.dotplot(result, var_names=result.var_names, groupby=["leiden"], save="leiden.png", dendrogram=True)
    # sc.pl.heatmap(result, var_names=result.var_names, groupby="leiden", save="leiden.png")
    return result
