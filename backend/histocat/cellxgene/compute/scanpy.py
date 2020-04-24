import scanpy as sc


def scanpy_umap(adata, obs_mask=None, pca_options={}, neighbors_options={}, umap_options={}):
    """
    Given adata and an obs mask, return a new embedding for adata[obs_mask, :]
    as an ndarray of shape (len(obs_mask), N), where N>=2.

    Do NOT mutate adata.
    """

    # backed mode is incompatible with the current implementation
    if adata.isbacked:
        raise NotImplementedError("Backed mode is incompatible with re-embedding")

    # https://github.com/theislab/anndata/issues/311
    obs_mask = slice(None) if obs_mask is None else obs_mask
    adata = adata[obs_mask, :].copy()

    for k in list(adata.obsm.keys()):
        del adata.obsm[k]
    for k in list(adata.uns.keys()):
        del adata.uns[k]

    sc.pp.pca(adata, zero_center=None, n_comps=min(adata.n_obs - 1, 50), **pca_options)
    sc.pp.neighbors(adata, **neighbors_options)
    sc.tl.umap(adata, **umap_options)

    return adata.obsm["X_umap"]
