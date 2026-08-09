"""Expression-source selection and projected extraction for plotting data.

annplyr is ggann's authoritative AnnData-to-tabular adapter.  Its public v0.3
accessor projects requested matrix columns before reading dense, sparse, view,
or backed inputs.  This module keeps source selection and ggann's user-facing
errors in one place while delegating every tabular matrix read to that accessor.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence

import annplyr as _annplyr  # registers the AnnData ``.ap`` accessor
import pandas as pd
from anndata import AnnData
from scipy import sparse


def resolve_source(
    adata,
    layer: str | None,
    use_raw: bool | None,
) -> tuple[str, str | None]:
    """Return the selected expression source as ``(kind, layer)``."""
    if layer is not None:
        if use_raw is True:
            raise ValueError("Cannot specify use_raw=True and a layer at the same time.")
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


def validate_materialization_budget(
    max_matrix_values: int | None,
    projected_values: int,
    *,
    context: str,
) -> None:
    """Reject an invalid or over-budget projected read before extraction.

    annplyr validates cumulative budgets within one accessor call.  A ggann
    aesthetic can span several layers and embeddings, which necessarily uses
    several accessor calls.  This preflight applies the same logical-cell
    accounting across the complete ggann request before the first call.
    """
    if max_matrix_values is None:
        return
    if max_matrix_values < 0:
        raise _annplyr.AnnplyrError("max_matrix_values must be non-negative or None")
    if projected_values > max_matrix_values:
        raise _annplyr.AnnplyrError(
            f"{context} would materialize {projected_values} matrix values, "
            f"which exceeds max_matrix_values={max_matrix_values}"
        )


def _validate_genes(
    adata,
    genes: Iterable[str],
    *,
    kind: str,
    layer: str | None,
) -> list[str]:
    """Validate and deduplicate a projected feature selection."""
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
    return genes


def _expression_kwargs(kind: str, layer: str | None, request) -> dict:
    """Map a resolved ggann source to public annplyr accessor arguments."""
    if kind == "raw":
        return {"raw": request}
    return {"x": request, "layer": layer if kind == "layer" else None}


def expression_frame(
    adata,
    genes: Iterable[str],
    *,
    kind: str,
    layer: str | None,
    max_matrix_values: int | None = None,
) -> tuple[pd.DataFrame, list[str]]:
    """Read only requested expression columns through public annplyr APIs.

    The returned frame retains pandas sparse columns when annplyr produces
    them.  Column names and rows follow the requested feature order and the
    input observation order, including duplicate observation names.
    """
    genes = _validate_genes(adata, genes, kind=kind, layer=layer)
    kwargs = _expression_kwargs(kind, layer, genes)
    frame = adata.ap.to_df(
        **kwargs,
        max_matrix_values=max_matrix_values,
    )
    # ``to_df(raw=...)`` prefixes columns to disambiguate mixed X/raw exports.
    # This single-source boundary is already unambiguous, so restore the public
    # feature names positionally rather than parsing names such as ``raw_raw_x``.
    frame.columns = genes
    return frame, genes


def expression_with_obs_frame(
    adata,
    genes: Iterable[str],
    obs: Sequence[str],
    *,
    kind: str,
    layer: str | None,
    max_matrix_values: int | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, list[str]]:
    """Project expression and grouping metadata in one public annplyr call.

    The two returned frames are split by position, so duplicate or reordered
    observation names remain safe.  A rare exported-name collision falls back
    to separate public calls, preserving the pre-existing behavior rather than
    exposing annplyr's combined-frame name-repair error.
    """
    genes = _validate_genes(adata, genes, kind=kind, layer=layer)
    obs = tuple(ordered_unique(obs))
    exported_genes = [f"raw_{gene}" for gene in genes] if kind == "raw" else genes
    exported_names = [*obs, *exported_genes]
    can_fuse = len(exported_names) == len(set(exported_names))

    if can_fuse:
        frame = adata.ap.to_df(
            obs=list(obs) if obs else None,
            **_expression_kwargs(kind, layer, genes),
            max_matrix_values=max_matrix_values,
        )
        split = len(obs)
        observations = frame.iloc[:, :split].copy(deep=False)
        expression = frame.iloc[:, split:].copy(deep=False)
        expression.columns = genes
        return expression, observations, genes

    expression, genes = expression_frame(
        adata,
        genes,
        kind=kind,
        layer=layer,
        max_matrix_values=max_matrix_values,
    )
    observations = adata.ap.to_df(obs=list(obs))
    return expression, observations, genes


def project_expression(
    adata,
    genes: Iterable[str],
    *,
    kind: str,
    layer: str | None,
    obs: Sequence[str] = (),
    max_matrix_values: int | None = None,
) -> tuple[AnnData, list[str]]:
    """Build a bounded temporary AnnData using public annplyr extraction.

    This private compatibility seam remains for extensions that need an
    AnnData-shaped projection.  New ggann code should prefer
    :func:`expression_frame` so it does not construct a temporary container.
    The input is never modified.
    """
    expression, genes = expression_frame(
        adata,
        genes,
        kind=kind,
        layer=layer,
        max_matrix_values=max_matrix_values,
    )

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

    sparse_columns = [isinstance(dtype, pd.SparseDtype) for dtype in expression.dtypes]
    matrix = (
        sparse.csr_matrix(expression.sparse.to_coo())
        if sparse_columns and all(sparse_columns)
        else expression.to_numpy(copy=False)
    )
    names = source_var_names(adata, kind)
    var = pd.DataFrame(index=pd.Index(genes, name=names.name))
    projected = AnnData(X=matrix, obs=obs_frame, var=var)
    return projected, genes


def densify_frame(frame: pd.DataFrame) -> pd.DataFrame:
    """Convert pandas sparse columns to dense columns in one block."""
    sparse_columns = [
        column for column in frame.columns if isinstance(frame[column].dtype, pd.SparseDtype)
    ]
    if not sparse_columns:
        return frame

    # pandas 3 represents scipy structural zeros with ``fill_value=NaN`` in
    # SparseArray while its sparse-to-COO conversion retains the real zeros.
    # Round-trip that block through scipy, then restore each column's numeric
    # subtype so float32 inputs do not silently widen.
    values = frame[sparse_columns].sparse.to_coo().toarray()
    dense_sparse = pd.DataFrame(
        {
            column: values[:, position].astype(
                frame[column].dtype.subtype,
                copy=False,
            )
            for position, column in enumerate(sparse_columns)
        },
        index=frame.index,
    )
    dense_columns = [column for column in frame.columns if column not in sparse_columns]
    if not dense_columns:
        return dense_sparse

    dense = pd.concat([frame[dense_columns], dense_sparse], axis=1)
    return dense.loc[:, frame.columns]
