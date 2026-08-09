"""Small deterministic AnnData fixtures for the executable vignettes."""

from __future__ import annotations

import hashlib

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import sparse

GENES = ["CD3D", "NKG7", "MS4A1", "CST3", "MKI67", "GAPDH"]
CELL_TYPES = ["T cell", "NK cell", "B cell", "Monocyte"]


def _matrix(values: np.ndarray, storage: str):
    if storage == "dense":
        return values.copy()
    if storage == "csr":
        return sparse.csr_matrix(values)
    if storage == "csc":
        return sparse.csc_matrix(values)
    raise ValueError(f"storage must be 'dense', 'csr', or 'csc', got {storage!r}")


def make_adata(*, storage: str = "dense", n_obs: int = 72) -> AnnData:
    """Return a reproducible single-cell fixture without downloading data."""
    rng = np.random.default_rng(7)
    group_index = np.arange(n_obs) % len(CELL_TYPES)
    counts = rng.poisson(0.8, size=(n_obs, len(GENES))).astype(np.float32)
    counts[np.arange(n_obs), group_index] += 5
    counts[np.arange(n_obs) % 2 == 1, 4] += 2
    expression = np.log1p(counts).astype(np.float32)

    obs = pd.DataFrame(
        {
            "cell_type": pd.Categorical(
                np.asarray(CELL_TYPES, dtype=object)[group_index],
                categories=CELL_TYPES,
                ordered=True,
            ),
            "condition": pd.Categorical(
                np.where(np.arange(n_obs) % 2 == 0, "control", "stimulated"),
                categories=["control", "stimulated"],
                ordered=True,
            ),
            "phase": pd.Categorical(
                np.asarray(["G1", "S", "G2M"], dtype=object)[np.arange(n_obs) % 3],
                categories=["G1", "S", "G2M"],
                ordered=True,
            ),
            "quality_score": rng.uniform(0.6, 1.0, size=n_obs),
        },
        index=pd.Index([f"cell_{i:03d}" for i in range(n_obs)], name="cell"),
    )
    var = pd.DataFrame(index=pd.Index(GENES, name="gene"))
    adata = AnnData(X=_matrix(expression, storage), obs=obs, var=var)
    adata.layers["counts"] = _matrix(counts, storage)
    adata.layers["logcounts"] = _matrix(expression, storage)
    adata.raw = AnnData(X=_matrix(expression, storage), var=var.copy())

    centers = np.asarray([[-2, 0], [0, 2], [2, 0], [0, -2]], dtype=np.float32)
    adata.obsm["X_umap"] = centers[group_index] + rng.normal(scale=0.28, size=(n_obs, 2)).astype(
        np.float32
    )
    adata.uns["cell_type_colors"] = ["#4477AA", "#EE6677", "#228833", "#CCBB44"]
    return adata


def fingerprint(adata: AnnData) -> str:
    """Hash stable public state used by the vignettes' mutation checks."""
    digest = hashlib.sha256()

    def update_matrix(matrix) -> None:
        if sparse.issparse(matrix):
            matrix = matrix.tocsr()
            for values in (matrix.data, matrix.indices, matrix.indptr):
                digest.update(np.ascontiguousarray(values).tobytes())
        else:
            digest.update(np.ascontiguousarray(np.asarray(matrix)).tobytes())

    update_matrix(adata.X)
    for key in sorted(key for key in adata.layers if key is not None):
        digest.update(str(key).encode())
        update_matrix(adata.layers[key])
    if adata.raw is not None:
        update_matrix(adata.raw.X)
    for key in sorted(adata.obsm):
        digest.update(key.encode())
        update_matrix(adata.obsm[key])
    digest.update("\0".join(map(str, adata.obs_names)).encode())
    digest.update("\0".join(map(str, adata.var_names)).encode())
    digest.update(adata.obs.astype(str).to_csv().encode())
    digest.update(repr(adata.uns).encode())
    return digest.hexdigest()
