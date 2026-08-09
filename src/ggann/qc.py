"""Quality-control plots over ``obs`` metrics and expression.

Observation metrics use the ``adata.ap`` extraction path. Highest-expression
ranking reads the resolved matrix directly so sparse and backed data remain
bounded during the whole-matrix calculation.
"""

from __future__ import annotations

import warnings
from typing import Sequence

import numpy as np
import pandas as pd
import plotnine_extra as pe
from plotnine import (
    aes,
    element_blank,
    geom_boxplot,
    geom_line,
    geom_point,
    geom_violin,
    ggplot,
    labs,
    scale_y_log10,
    theme,
)

from ._expression import densify_frame as _densify
from ._expression import resolve_source, source_label, source_matrix
from ._palette import scale_fill_obs
from .theme import theme_ggann

__all__ = [
    "plot_qc_violin",
    "plot_qc_scatter",
    "plot_highest_expr_genes",
    "plot_variance_ratio",
]

# Common QC metric names across scanpy / older workflows.
_DEFAULT_METRICS = [
    "n_genes_by_counts",
    "n_genes",
    "total_counts",
    "n_counts",
    "pct_counts_mt",
    "pct_counts_ribo",
    "percent_mito",
]


def _resolve_metrics(adata, metrics):
    if metrics is not None:
        missing = [m for m in metrics if m not in adata.obs]
        if missing:
            raise KeyError(f"QC metrics not in adata.obs: {missing}")
        return list(metrics)
    found = [m for m in _DEFAULT_METRICS if m in adata.obs]
    if not found:
        raise ValueError("No default QC metrics found in adata.obs; pass metrics=[...] explicitly.")
    return found


def _dense_row_blocks(matrix, *, target_elements: int = 1_000_000):
    """Yield bounded dense row blocks from an in-memory or backed matrix."""
    rows_per_block = max(1, target_elements // max(matrix.shape[1], 1))
    for start in range(0, matrix.shape[0], rows_per_block):
        stop = min(start + rows_per_block, matrix.shape[0])
        yield start, stop, np.asarray(matrix[start:stop])


def _streamed_dense_totals(matrix) -> tuple[np.ndarray, bool]:
    """Compute NaN-skipping totals and negativity without loading a backed matrix."""
    totals = np.empty(matrix.shape[0], dtype=float)
    has_negative = False
    with np.errstate(invalid="ignore"):
        for start, stop, values in _dense_row_blocks(matrix):
            totals[start:stop] = np.nansum(values, axis=1)
            has_negative = has_negative or bool(np.any(values < 0))
    return totals, has_negative


def _expression_totals(matrix) -> tuple[np.ndarray, bool]:
    """Return row sums that skip NaNs and whether element-wise means are needed."""
    from scipy import sparse

    totals = np.asarray(matrix.sum(axis=1)).reshape(-1).astype(float, copy=False)
    has_nan_sum = bool(np.isnan(totals).any())
    if has_nan_sum:
        if sparse.issparse(matrix):
            cleaned = matrix.copy()
            if np.issubdtype(cleaned.data.dtype, np.inexact):
                cleaned.data[np.isnan(cleaned.data)] = 0
            with np.errstate(invalid="ignore"):
                totals = np.asarray(cleaned.sum(axis=1)).reshape(-1).astype(float, copy=False)
        else:
            totals = np.empty(matrix.shape[0], dtype=float)
            with np.errstate(invalid="ignore"):
                for start, stop, values in _dense_row_blocks(matrix):
                    totals[start:stop] = np.nansum(values, axis=1)
    needs_elementwise = has_nan_sum or not bool(np.isfinite(totals).all())
    return totals, needs_elementwise


def _mean_percentages(
    matrix,
    totals: np.ndarray,
    valid_totals: np.ndarray,
    *,
    elementwise: bool,
) -> np.ndarray:
    """Mean per-gene percentages, excluding undefined values like pandas."""
    from scipy import sparse

    n_vars = matrix.shape[1]
    if not elementwise:
        if not valid_totals.any():
            return np.full(n_vars, np.nan, dtype=float)
        weights = np.zeros_like(totals, dtype=float)
        weights[valid_totals] = 1.0 / totals[valid_totals]
        return np.asarray(matrix.T @ weights).reshape(-1) * (100.0 / valid_totals.sum())

    means = np.full(n_vars, np.nan, dtype=float)
    if sparse.issparse(matrix):
        columns = matrix.tocsc(copy=True)
        columns.sum_duplicates()
        denominator = int(valid_totals.sum())
        for index in range(n_vars):
            start, stop = columns.indptr[index : index + 2]
            rows = columns.indices[start:stop]
            keep = valid_totals[rows]
            fractions = np.full(stop - start, np.nan, dtype=float)
            with np.errstate(divide="ignore", invalid="ignore"):
                np.divide(
                    columns.data[start:stop],
                    totals[rows],
                    out=fractions,
                    where=keep,
                )
            defined = keep & ~np.isnan(fractions)
            count = denominator - int((keep & ~defined).sum())
            if count:
                with np.errstate(invalid="ignore"):
                    means[index] = float(fractions[defined].sum()) * (100.0 / count)
        return means

    numerators = np.zeros(n_vars, dtype=float)
    counts = np.zeros(n_vars, dtype=np.int64)
    for start, stop, values in _dense_row_blocks(matrix):
        fractions = np.full(values.shape, np.nan, dtype=float)
        with np.errstate(divide="ignore", invalid="ignore"):
            np.divide(
                values,
                totals[start:stop, None],
                out=fractions,
                where=valid_totals[start:stop, None],
            )
        defined = ~np.isnan(fractions)
        with np.errstate(invalid="ignore"):
            numerators += np.nansum(fractions, axis=0)
        counts += defined.sum(axis=0)
    np.divide(numerators, counts, out=means, where=counts > 0)
    means *= 100.0
    return means


def _select_backed_dense_columns(matrix, indices: np.ndarray) -> np.ndarray:
    """Read selected backed columns in bounded row blocks and requested order."""
    result = np.empty((matrix.shape[0], len(indices)), dtype=matrix.dtype)
    if not len(indices):
        return result
    order = np.argsort(indices)
    inverse = np.argsort(order)
    sorted_indices = indices[order]
    rows_per_block = max(1, 1_000_000 // len(indices))
    for start in range(0, matrix.shape[0], rows_per_block):
        stop = min(start + rows_per_block, matrix.shape[0])
        values = np.asarray(matrix[start:stop, sorted_indices])
        result[start:stop] = values[:, inverse]
    return result


def plot_qc_violin(
    adata,
    metrics: Sequence[str] | None = None,
    group_by: str | None = None,
    *,
    scale: str = "width",
):
    """Violin distributions of per-cell QC metrics, one facet per metric.

    With ``group_by`` the violins are split by that obs column (stored palette
    reused); without it, one violin per metric.

    Parameters
    ----------
    adata
        Annotated data matrix.
    metrics : sequence of str, optional
        Numeric observation columns; common QC fields are detected by default.
    group_by : str, optional
        Observation column defining violins.
    scale : str
        plotnine violin-width scaling mode.

    Returns
    -------
    plotnine.ggplot
        Composable QC violin plot.

    Raises
    ------
    KeyError
        If an explicit metric or grouping column is absent.
    ValueError
        If no default QC metric exists or scaling is invalid.

    Notes
    -----
    Only observation metadata is copied into the plot; ``adata`` is unchanged.

    Examples
    --------
    >>> p = plot_qc_violin(adata, metrics=["total_counts", "pct_counts_mt"])
    """
    metrics = _resolve_metrics(adata, metrics)
    cols = ([group_by] if group_by else []) + metrics
    df = _densify(adata.ap.to_df(obs=cols))
    id_vars = [group_by] if group_by else []
    long = df.melt(id_vars=id_vars, value_vars=metrics, var_name="metric", value_name="value")

    if group_by:
        p = (
            ggplot(long, aes(group_by, "value", fill=group_by))
            + geom_violin(scale=scale)
            + scale_fill_obs(adata, group_by)
        )
    else:
        long["_x"] = "all cells"
        p = ggplot(long, aes("_x", "value")) + geom_violin(scale=scale, fill="#4c72b0")

    p = (
        p
        + pe.facet_wrap("~metric", scales="free_y")
        + labs(x="", y="", fill=group_by or "")
        + theme_ggann()
    )
    if group_by:
        # the fill legend already encodes group; many long cell-type names across
        # shared-x facets are illegible -- drop the redundant x ticks (added last so
        # it wins over theme_ggann's axis_text).
        p = p + theme(axis_text_x=element_blank(), axis_ticks_major_x=element_blank())
    return p


def plot_qc_scatter(
    adata,
    x: str,
    y: str,
    color: str | None = None,
    *,
    layer: str | None = None,
    use_raw: bool | None = None,
    size: float = 1.0,
):
    """Scatter of two obs QC metrics (e.g. total_counts vs pct_counts_mt).

    Thin wrapper over the grammar path; ``color`` may be an obs column or a gene.

    Parameters
    ----------
    adata
        Annotated data matrix.
    x, y : str or selector
        Observation columns or expression variables.
    color : str or selector, optional
        Observation column or expression variable.
    layer, use_raw : optional
        Mutually exclusive expression source.
    size : float
        Point size.

    Returns
    -------
    plotnine.ggplot
        Composable QC scatter plot.

    Raises
    ------
    KeyError
        If an explicit field or layer is missing.
    ValueError
        If expression-source selection is invalid.

    Notes
    -----
    Only mapped fields are projected through the grammar path; input is unchanged.

    Examples
    --------
    >>> p = plot_qc_scatter(adata, x="total_counts", y="n_genes_by_counts")
    """
    from plotnine import geom_point

    from ._palette import scale_color_obs
    from ._resolve import plain_name
    from .grammar import aes as _aes
    from .grammar import gganndata

    mapping = _aes(x, y) if color is None else _aes(x, y, color=color)
    plot = gganndata(adata, mapping, layer=layer, use_raw=use_raw)
    fields = [x, y, *([color] if color is not None else [])]
    missing = [field for field in fields if plain_name(adata, field) not in plot.data.columns]
    if missing:
        raise KeyError(f"Could not resolve field(s) from observations or expression: {missing}.")

    plot = plot + geom_point(size=size, alpha=0.6)
    if color is not None:
        color_name = plain_name(adata, color)
        is_categorical_obs = color_name in adata.obs and isinstance(
            adata.obs[color_name].dtype, pd.CategoricalDtype
        )
        if is_categorical_obs:
            plot = plot + scale_color_obs(adata, color_name)
        elif pd.api.types.is_numeric_dtype(plot.data[color_name]):
            # numeric (gene or continuous metric) -> expression colourmap, matching plot_embedding
            from .theme import scale_color_expression

            plot = plot + scale_color_expression()
    return plot


def plot_highest_expr_genes(adata, n: int = 20, *, use_raw: bool = False, layer: str | None = None):
    """Boxplots of the ``n`` genes accounting for the most counts per cell.

    Like ``sc.pl.highest_expr_genes``, this reads ``adata.X`` by default (pass
    ``layer="counts"`` or ``use_raw=True`` to point elsewhere): each cell's values
    are turned into percentages of that cell's total, and genes are ranked by mean
    per-cell percentage. For meaningful results ``adata.X`` should hold counts or
    normalized (not scaled) expression. This whole-matrix calculation stays sparse
    until the selected genes are prepared for plotting.

    Parameters
    ----------
    adata
        Annotated data matrix.
    n : int
        Number of top genes.
    use_raw : bool
        Read ``adata.raw`` rather than ``adata.X``.
    layer : str, optional
        Named expression layer.

    Returns
    -------
    plotnine.ggplot
        Composable percentage box plot.

    Raises
    ------
    KeyError
        If the selected layer is absent.
    ValueError
        If ``n`` is non-positive or source selection is invalid.

    Notes
    -----
    Ranking requires a streamed or sparse pass over the selected whole matrix by
    definition. Only top-gene columns are then materialized; input is unchanged.

    Examples
    --------
    >>> p = plot_highest_expr_genes(adata, n=10, layer="counts")
    """
    from scipy import sparse

    if isinstance(n, (bool, np.bool_)) or not isinstance(n, (int, np.integer)) or n < 1:
        raise ValueError(f"n must be a positive integer, got {n!r}.")
    kind, lyr = resolve_source(adata, layer, use_raw)
    matrix, var_names = source_matrix(adata, kind, lyr)
    if not var_names.is_unique:
        raise ValueError(
            f"Variable names in {source_label(kind, lyr)} must be unique before plotting."
        )
    if hasattr(matrix, "to_memory"):
        matrix = matrix.to_memory()

    streamed_dense = not sparse.issparse(matrix) and not hasattr(matrix, "sum")
    if streamed_dense:
        totals, has_negative = _streamed_dense_totals(matrix)
        elementwise = True
    else:
        has_negative = (
            bool(matrix.data.size and np.any(matrix.data < 0))
            if sparse.issparse(matrix)
            else bool(np.any(np.asarray(matrix) < 0))
        )
        totals, elementwise = _expression_totals(matrix)
    if has_negative:
        warnings.warn(
            "plot_highest_expr_genes: the expression matrix has negative values, which "
            "looks like scaled/z-scored data -- 'percent of total counts' will be "
            "meaningless. Pass use_raw=True or layer= to point at counts or "
            "log-normalized values.",
            stacklevel=2,
        )
    zero_total = totals == 0
    valid = ~np.isnan(totals) & ~zero_total
    n_zero = int(zero_total.sum())
    if n_zero:
        warnings.warn(
            f"{n_zero} cell(s) with zero total counts excluded from plot_highest_expr_genes.",
            stacklevel=2,
        )
    means = _mean_percentages(
        matrix,
        totals,
        valid,
        elementwise=elementwise,
    )
    n = min(n, len(var_names))
    top_indices = np.argsort(-means, kind="stable")[:n]
    top = var_names[top_indices].tolist()

    if streamed_dense:
        selected = _select_backed_dense_columns(matrix, top_indices)
    else:
        selected = matrix[:, top_indices]
        selected = selected.toarray() if sparse.issparse(selected) else np.asarray(selected)
    percentages = np.full(selected.shape, np.nan, dtype=float)
    np.divide(
        selected,
        totals[:, None],
        out=percentages,
        where=valid[:, None],
    )
    percentages *= 100.0
    long = pd.DataFrame(percentages, columns=top).melt(var_name="gene", value_name="percent")
    long["gene"] = pd.Categorical(long["gene"], categories=list(reversed(top)), ordered=True)
    return (
        ggplot(long, aes("gene", "percent"))
        + geom_boxplot(fill="#4c72b0", outlier_alpha=0.2)
        + pe.coord_flip()
        + labs(x="", y="% of total counts per cell")
        + theme_ggann()
    )


def plot_variance_ratio(adata, n_pcs: int = 50, *, key: str = "pca", log: bool = True):
    """Scree / elbow plot of PCA variance ratio (``sc.pl.pca_variance_ratio``).

    Reads the per-PC variance ratio stored by ``sc.tl.pca`` in
    ``adata.uns[key]['variance_ratio']`` and draws it as a point-and-line elbow, to
    help pick how many PCs to keep. ``log=True`` (scanpy's default) puts the y axis
    on a log scale so the elbow is easier to read.

    Parameters
    ----------
    adata
        Annotated data matrix with stored PCA results.
    n_pcs : int
        Maximum principal components to display.
    key : str
        Key in ``adata.uns``.
    log : bool
        Use a logarithmic y scale.

    Returns
    -------
    plotnine.ggplot
        Composable scree plot.

    Raises
    ------
    KeyError
        If variance-ratio information is absent.
    ValueError
        If ``n_pcs`` is invalid.

    Notes
    -----
    Reads only a small vector from ``adata.uns`` and does not mutate input.

    Examples
    --------
    >>> p = plot_variance_ratio(adata, n_pcs=20)
    """
    uns = adata.uns.get(key)
    if uns is None or "variance_ratio" not in uns:
        raise KeyError(
            f"adata.uns['{key}']['variance_ratio'] not found; run sc.tl.pca(adata) first."
        )
    vr = np.asarray(uns["variance_ratio"], dtype=float)
    n = min(n_pcs, vr.size)
    df = pd.DataFrame({"PC": np.arange(1, n + 1), "variance_ratio": vr[:n]})
    plot = (
        ggplot(df, aes("PC", "variance_ratio"))
        + geom_line(color="#b0b0b0")
        + geom_point(color="#2166ac", size=1.8)
        + labs(x="principal component", y="variance ratio")
        + theme_ggann()
    )
    if log:
        plot = plot + scale_y_log10()
    return plot
