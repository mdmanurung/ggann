"""Shared PBMC3K loader and state fingerprint for the executable vignettes.

Every vignette runs on the same real dataset: ``scanpy.datasets.pbmc3k_processed``,
the 2,638-cell 10x Genomics PBMC sample used by Scanpy's and Seurat's clustering
tutorials. Readers who have followed either tutorial recognise its Louvain cell
types and marker genes immediately.

The published object is left intact. :func:`load_adata` only adds layers and two
documented recodings of columns the dataset already carries, so ``X_umap``,
``X_pca``, and ``rank_genes_groups`` keep corresponding to the matrix they were
computed from.
"""

from __future__ import annotations

import hashlib
import os
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData
from scipy import sparse

#: Cache location for downloaded Scanpy datasets. Pinned to the repository so a
#: vignette resolves the same file regardless of the working directory.
DATASET_DIR = Path(__file__).resolve().parents[2] / "data"

#: Filename Scanpy writes for ``pbmc3k_processed`` inside ``DATASET_DIR``.
DATASET_FILE = "pbmc3k_processed.h5ad"

#: Louvain labels published with the dataset, ordered lymphoid before myeloid.
CELL_TYPES = [
    "CD4 T cells",
    "CD8 T cells",
    "NK cells",
    "B cells",
    "CD14+ Monocytes",
    "FCGR3A+ Monocytes",
    "Dendritic cells",
    "Megakaryocytes",
]

#: Lineage compartments assigned to those labels by standard haematopoiesis.
LYMPHOID = ("CD4 T cells", "CD8 T cells", "NK cells", "B cells")

#: Recognisable PBMC markers that survived the dataset's highly-variable-gene
#: selection, so they are readable from ``.X`` without falling back to ``.raw``.
#: ``CD3D``, ``IL7R``, and ``LYZ`` did not survive it; reach them with
#: ``use_raw=True``.
MARKERS = ["NKG7", "GNLY", "MS4A1", "CST3"]


def _matrix(values, storage: str):
    if storage == "dense":
        return np.asarray(values.todense()) if sparse.issparse(values) else np.asarray(values)
    if storage == "csr":
        return sparse.csr_matrix(values)
    if storage == "csc":
        return sparse.csc_matrix(values)
    raise ValueError(f"storage must be 'dense', 'csr', or 'csc', got {storage!r}")


def dataset_path() -> Path:
    """Return the cached dataset path, or explain how to populate it."""
    path = DATASET_DIR / DATASET_FILE
    if os.environ.get("GGANN_DOCS_OFFLINE") == "1" and not path.exists():
        raise FileNotFoundError(
            f"{path} is missing and GGANN_DOCS_OFFLINE=1 forbids downloading it. "
            "Run `python scripts/fetch_datasets.py` once while online."
        )
    return path


def load_adata(*, storage: str = "csr") -> AnnData:
    """Return the processed PBMC3K dataset prepared for the vignettes.

    Parameters
    ----------
    storage
        Container for ``.X`` and ``layers['logcounts']``: ``'dense'``, ``'csr'``,
        or ``'csc'``. ``layers['scaled']`` stays dense because a z-scored matrix
        has no zeros worth storing sparsely.

    Notes
    -----
    Three additions to the published object, all derived from columns it already
    carries:

    ``layers['scaled']``
        The published ``.X``: highly-variable genes, regressed and z-scored. This
        is the matrix behind ``X_pca``, ``X_umap``, and ``rank_genes_groups``.
    ``layers['logcounts']`` and ``.X``
        Log-normalised expression for the same genes, read from ``.raw``. Marker
        plots want a non-negative scale, so this becomes the default matrix.
    ``obs['compartment']`` and ``obs['depth']``
        A lymphoid/myeloid recoding of ``obs['louvain']`` and a median split of
        the measured ``obs['n_counts']``.

    Two presentation changes go beyond adding fields. ``obs['louvain']`` is
    reordered to :data:`CELL_TYPES`, which groups the lymphoid labels first, and
    ``uns['louvain_colors']`` is replaced with a palette matching that order.
    ``uns['rank_genes_groups']`` still carries the published category order, so
    a workflow that reads it must realign rather than assume positional
    correspondence.
    """
    dataset_path()
    sc.settings.datasetdir = DATASET_DIR
    adata = sc.datasets.pbmc3k_processed()

    adata.layers["scaled"] = np.asarray(adata.X)
    lognormalised = adata.raw[:, adata.var_names].X
    adata.layers["logcounts"] = _matrix(lognormalised, storage)
    adata.X = _matrix(lognormalised, storage)

    adata.obs["louvain"] = pd.Categorical(
        adata.obs["louvain"].astype(str), categories=CELL_TYPES, ordered=True
    )
    adata.obs["compartment"] = pd.Categorical(
        np.where(adata.obs["louvain"].isin(LYMPHOID), "Lymphoid", "Myeloid"),
        categories=["Lymphoid", "Myeloid"],
        ordered=True,
    )
    adata.obs["depth"] = pd.Categorical(
        np.where(adata.obs["n_counts"] < adata.obs["n_counts"].median(), "low", "high"),
        categories=["low", "high"],
        ordered=True,
    )
    adata.uns["louvain_colors"] = np.asarray(
        ["#4477AA", "#66CCEE", "#228833", "#CCBB44", "#EE6677", "#AA3377", "#BBBBBB", "#332288"]
    )
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
    # Scanpy plotting caches palettes under ``uns['<key>_colors']``; hashing the
    # keys and those palettes is what makes the mutation checks meaningful.
    for key in sorted(adata.uns):
        digest.update(key.encode())
        if key.endswith("_colors"):
            digest.update(np.asarray(adata.uns[key]).astype(str).tobytes())
    return digest.hexdigest()
