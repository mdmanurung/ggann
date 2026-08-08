"""Expression-source selection and projection for plotting data.

annplyr is ggann's tabular extraction adapter.  Its selectors operate after a
matrix has been represented as a pandas frame, so passing a full sparse
``AnnData`` makes the cost scale with every variable rather than the requested
genes.  This module projects the matrix first and then hands the small,
aligned object to annplyr.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence

import annplyr as _annplyr  # noqa: F401 -- registers the AnnData ``.ap`` accessor
from anndata import AnnData
import numpy as np
import pandas as pd


def resolve_source(
    adata,
    layer: str | None,
    use_raw: bool | None,
) -> tuple[str, str | None]:
    """Return the selected expression source as ``(kind, layer)``."""
    if layer is not None:
        if use_raw is True:
            raise ValueError(
                "Cannot specify use_raw=True and a layer at the same time."
            )
        if layer not in adata.layers:
            raise KeyError(f"Layer {layer!r} not found in adata.layers.")
        return "layer", layer
    if use_raw is None:
        use_raw = adata.raw is not None
    if use_raw:
        if adata.raw is None:
            raise ValueError("use_raw=True but adata.raw is None.")
        return "raw", None
    return "x", None


def source_var_names(adata, kind: str) -> pd.Index:
    """Return variable names for a resolved source without copying them."""
    return adata.raw.var_names if kind == "raw" else adata.var_names


def source_label(kind: str, layer: str | None) -> str:
    """Return a user-facing expression-source label."""
    return {"raw": ".raw", "x": ".X"}.get(kind, f"layer {layer!r}")


def source_matrix(adata, kind: str, layer: str | None):
    """Return a resolved matrix and its variable names.

    Whole-matrix algorithms use this private seam so direct AnnData container
    access does not leak into plotting modules.
    """
    if kind == "raw":
        return adata.raw.X, adata.raw.var_names
    if kind == "layer":
        return adata.layers[layer], adata.var_names
    return adata.X, adata.var_names


def ordered_unique(values: Iterable[str]) -> list[str]:
    """Return values in first-seen order, with duplicates removed."""
    return list(dict.fromkeys(values))


def project_expression(
    adata,
    genes: Iterable[str],
    *,
    kind: str,
    layer: str | None,
    obs: Sequence[str] = (),
) -> tuple[AnnData, list[str]]:
    """Project one expression source to requested genes and observation data.

    The returned ``AnnData`` contains only the matrix columns and observation
    columns needed by the caller.  Matrix selection happens before annplyr
    builds a pandas representation, avoiding work proportional to unrelated
    variables.  The input object is never modified.
    """
    genes = ordered_unique(genes)
    names = source_var_names(adata, kind)
    if not names.is_unique:
        raise ValueError(
            f"Variable names in {source_label(kind, layer)} must be unique before plotting."
        )

    missing = [gene for gene in genes if gene not in names]
    if missing:
        quoted = ", ".join(repr(gene) for gene in missing)
        raise KeyError(f"Gene(s) not found in {source_label(kind, layer)}: {quoted}.")

    obs = tuple(ordered_unique(obs))
    if obs:
        missing_obs = [name for name in obs if name not in adata.obs.columns]
        if missing_obs:
            quoted = ", ".join(repr(name) for name in missing_obs)
            raise KeyError(f"Observation column(s) not found: {quoted}.")
        obs_frame = adata.ap.to_df(obs=list(obs))
    else:
        obs_frame = pd.DataFrame(index=adata.obs_names.copy())

    # annplyr joins tidy expression and observation data by ``obs_name``. Duplicate
    # caller names would therefore create a many-to-many expansion in the temporary
    # projection. Give only that private object positional names; callers restore
    # the original labels after extraction.
    if not obs_frame.index.is_unique:
        obs_frame.index = pd.Index(
            [str(index) for index in range(len(obs_frame))],
            name=obs_frame.index.name,
        )

    if kind == "raw":
        # Raw backed views in anndata 0.12 expose the unsliced parent matrix.
        # Slice the backing matrix directly, with sorted indices for h5py, then
        # restore the caller's first-seen feature order on the bounded result.
        indices = names.get_indexer(genes)
        order = np.argsort(indices)
        matrix = adata.raw.X[:, indices[order]]
        if not np.array_equal(order, np.arange(len(order))):
            matrix = matrix[:, np.argsort(order)]
    else:
        view = adata[:, genes]
        matrix = view.layers[layer] if kind == "layer" else view.X

    var = pd.DataFrame(index=pd.Index(genes, name=names.name))
    projected = AnnData(X=matrix, obs=obs_frame, var=var)
    return projected, genes


def densify_frame(frame: pd.DataFrame) -> pd.DataFrame:
    """Convert pandas sparse columns to dense columns in one block."""
    sparse_columns = [
        column
        for column in frame.columns
        if isinstance(frame[column].dtype, pd.SparseDtype)
    ]
    if not sparse_columns:
        return frame

    dense_sparse = frame[sparse_columns].sparse.to_dense()
    dense_columns = [column for column in frame.columns if column not in sparse_columns]
    if not dense_columns:
        return dense_sparse

    dense = pd.concat([frame[dense_columns], dense_sparse], axis=1)
    return dense.loc[:, frame.columns]
