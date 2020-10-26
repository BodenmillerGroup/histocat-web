from typing import Any, Dict

import scanpy as sc
from anndata import AnnData


def process(adata: AnnData, step: Dict[str, Any], output: Dict[str, Any]):
    """Cluster cells using the Louvain algorithm"""
    output["louvain"] = True

    resolution = step.get("resolution")
    flavor = step.get("flavor")
    directed = step.get("directed")
    use_weights = step.get("useWeights")
    result = sc.tl.louvain(
        adata, resolution=resolution, flavor=flavor, directed=directed, use_weights=use_weights, copy=True
    )
    sc.pl.stacked_violin(result, var_names=result.var_names, groupby=["louvain"], save="louvain.png")
    # sc.pl.dotplot(result, var_names=result.var_names, groupby=["louvain"], save="louvain.png")
    # sc.pl.heatmap(result, var_names=result.var_names, groupby=["louvain"], save="louvain.png")
    return result
