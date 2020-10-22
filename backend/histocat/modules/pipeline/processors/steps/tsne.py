from typing import Any, Dict

from anndata import AnnData
import scanpy as sc


def process(adata: AnnData, step: Dict[str, Any], output: Dict[str, Any]):
    """Calculate t-Distributed Stochastic Neighbor Embedding data"""
    output["tsne"] = True

    perplexity = step.get("perplexity")
    early_exaggeration = step.get("earlyExaggeration")
    learning_rate = step.get("learningRate")
    result = sc.tl.tsne(adata, n_pcs=2, perplexity=perplexity, early_exaggeration=early_exaggeration, learning_rate=learning_rate, copy=True)
    return result
