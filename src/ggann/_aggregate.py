"""Grouped expression aggregation over annplyr-projected expression frames.

Produces the two quantities scanpy's dotplot / matrixplot are built from:

* ``mean_expression`` -- mean of the (raw / X / layer) matrix per group per gene
* ``fraction``        -- fraction of cells with expression above a cutoff

Matching scanpy, aggregation defaults to ``adata.raw`` when present. With
``standard_scale=None`` (the dotplot default) the mean is the raw group mean, so
the numbers reproduce ``sc.pl.dotplot`` / ``sc.pl.matrixplot`` (which also
default to ``standard_scale=None``).
"""

from __future__ import annotations

from typing import Iterable, Literal

import numpy as np
import pandas as pd
from scipy import sparse

from ._expression import densify_frame as _densify
from ._expression import expression_with_obs_frame, ordered_unique, resolve_source

__all__ = [
    "aggregate_expression",
    "aggregate_means",
    "expression_source",
    "group_means",
    "tidy_expression",
]

# Internal compatibility name used by older tests and extensions.
expression_source = resolve_source

StandardScale = Literal["var", "group", "zscore"]


def _projected_frames(adata, genes, by, *, kind, layer):
    """Extract one expression projection and its grouping columns positionally."""
    expression, groups, genes = expression_with_obs_frame(
        adata,
        genes,
        by,
        kind=kind,
        layer=layer,
    )
    if len(expression) != len(groups):  # defensive: both accessors must preserve rows
        raise RuntimeError(
            "annplyr returned expression and observation projections with different row counts."
        )
    return expression.reset_index(drop=True), groups.reset_index(drop=True), genes


def _matrix_from_frame(frame: pd.DataFrame):
    """Return a dense array or sparse CSR matrix without widening the projection."""
    sparse_columns = [isinstance(dtype, pd.SparseDtype) for dtype in frame.dtypes]
    if sparse_columns and all(sparse_columns):
        return frame.sparse.to_coo().tocsr()
    if any(sparse_columns):
        frame = _densify(frame)
    return frame.to_numpy(copy=False)


def _grouped_expression(
    adata,
    genes,
    by,
    *,
    kind,
    layer,
    expression_cutoff: float | None,
) -> tuple[pd.DataFrame, pd.DataFrame | None]:
    """Compute group means and optional detection fractions in one vectorized pass."""
    expression, groups, genes = _projected_frames(
        adata,
        genes,
        by,
        kind=kind,
        layer=layer,
    )

    # Grouped summaries omit rows missing any grouping key, matching pandas and
    # the documented plotting contract. The projected expression row order is
    # identical, so positional filtering remains safe with duplicate obs names.
    valid = groups[list(by)].notna().all(axis=1).to_numpy()
    valid_groups = groups.loc[valid, list(by)].reset_index(drop=True)
    keys = pd.MultiIndex.from_frame(valid_groups, names=list(by))
    codes, _ = pd.factorize(keys, sort=False)
    first = ~keys.duplicated()
    group_frame = valid_groups.loc[first].reset_index(drop=True)
    n_groups = len(group_frame)

    matrix = _matrix_from_frame(expression)
    matrix = matrix[np.flatnonzero(valid)]
    counts = np.bincount(codes, minlength=n_groups).astype(float, copy=False)
    indicator = sparse.csr_matrix(
        (np.ones(len(codes), dtype=float), (codes, np.arange(len(codes)))),
        shape=(n_groups, len(codes)),
    )

    def _dense_product(values) -> np.ndarray:
        product = indicator @ values
        if sparse.issparse(product):
            product = product.toarray()
        return np.asarray(product, dtype=float)

    if n_groups:
        if sparse.issparse(matrix):
            missing = matrix.copy()
            missing.data = np.isnan(missing.data).astype(np.float32, copy=False)
            missing.eliminate_zeros()
            valid_counts = counts[:, None] - _dense_product(missing)
            clean = matrix.copy()
            clean.data = np.where(np.isnan(clean.data), 0, clean.data)
        else:
            missing = np.isnan(matrix)
            valid_counts = _dense_product(~missing)
            clean = np.where(missing, 0, matrix)
        sums = _dense_product(clean)
        means = np.full(sums.shape, np.nan, dtype=float)
        np.divide(sums, valid_counts, out=means, where=valid_counts > 0)
    else:
        means = np.empty((0, len(genes)), dtype=float)
    group_index = (
        pd.Index(group_frame[by[0]], name=by[0])
        if len(by) == 1
        else pd.MultiIndex.from_frame(group_frame, names=list(by))
    )
    mean_df = pd.DataFrame(means, index=group_index, columns=genes)

    if expression_cutoff is None:
        return mean_df, None

    if sparse.issparse(matrix):
        compared = matrix.copy()
        if expression_cutoff >= 0:
            compared.data = (compared.data > expression_cutoff).astype(np.float32, copy=False)
            compared.eliminate_zeros()
            detected = _dense_product(compared)
        else:
            # Implicit sparse zeros exceed a negative cutoff. Count only explicit
            # values that fail (including NaN), then subtract their group rate
            # from one.
            compared.data = (np.isnan(compared.data) | (compared.data <= expression_cutoff)).astype(
                np.float32, copy=False
            )
            compared.eliminate_zeros()
            failed = _dense_product(compared)
            detected = counts[:, None] - failed
    else:
        detected = _dense_product(np.greater(matrix, expression_cutoff))
    fractions = detected / counts[:, None] if n_groups else np.empty((0, len(genes)), dtype=float)
    fraction_df = pd.DataFrame(fractions, index=group_index, columns=genes)
    return mean_df, fraction_df


def tidy_expression(adata, genes, group_by, *, layer=None, use_raw=None, extra_obs=()):
    """Long per-cell expression ``[obs_name, feature, value, group_by, *extra_obs]``.

    Shared by the violin / stacked-violin / tracksplot paths so the source
    picking, densification and feature-ordering live in one place. ``extra_obs``
    carries additional obs columns (e.g. a ``split_by`` facet) into the frame.
    """
    genes = ordered_unique(genes)
    obs = list(dict.fromkeys([group_by, *extra_obs]))  # dedupe, keep order
    kind, lyr = resolve_source(adata, layer, use_raw)
    expression, observations, genes = _projected_frames(
        adata,
        genes,
        obs,
        kind=kind,
        layer=lyr,
    )
    values = _densify(expression).to_numpy(copy=False)
    n_genes = len(genes)
    tidy = pd.DataFrame(
        {
            "obs_name": pd.Series(adata.obs_names).repeat(n_genes).reset_index(drop=True),
            "feature": pd.Categorical(
                np.tile(genes, adata.n_obs),
                categories=genes,
                ordered=True,
            ),
            "value": values.reshape(-1),
        }
    )
    for name in obs:
        tidy[name] = observations[name].repeat(n_genes).reset_index(drop=True)
    return tidy


def group_means(
    adata,
    genes: Iterable[str],
    group_by: str,
    *,
    layer=None,
    use_raw=None,
    extra_by: Iterable[str] = (),
) -> pd.DataFrame:
    """Mean expression per group (index) per gene (columns), mean-only.

    Like :func:`aggregate_expression` but skips the fraction-expressing pass, for
    callers (e.g. :func:`ggann.plot_correlation`) that only need the group means.
    """
    genes = ordered_unique(genes)
    kind, lyr = resolve_source(adata, layer, use_raw)
    by = list(dict.fromkeys([group_by, *extra_by]))
    mean_df, _ = _grouped_expression(
        adata,
        genes,
        by,
        kind=kind,
        layer=lyr,
        expression_cutoff=None,
    )
    return mean_df


def _standardize(mean_df: pd.DataFrame, standard_scale: str | None) -> pd.DataFrame:
    """Optionally rescale mean expression across groups/genes.

    ``'var'`` (per-gene 0..1) and ``'group'`` (per-group 0..1) match scanpy's
    ``standard_scale``. ``'zscore'`` (per-gene z-score, population std) is an
    ggann extension not present in scanpy.
    """
    if standard_scale is None:
        return mean_df
    if standard_scale == "var":  # per-gene (column) 0..1
        rng = (mean_df.max() - mean_df.min()).replace(0, 1)
        return (mean_df - mean_df.min()) / rng
    if standard_scale == "group":  # per-group (row) 0..1
        rng = (mean_df.max(axis=1) - mean_df.min(axis=1)).replace(0, 1)
        return mean_df.sub(mean_df.min(axis=1), axis=0).div(rng, axis=0)
    if standard_scale == "zscore":  # per-gene z-score
        std = mean_df.std(ddof=0).replace(0, 1)
        return (mean_df - mean_df.mean()) / std
    raise ValueError(
        f"standard_scale must be None, 'var', 'group' or 'zscore', got {standard_scale!r}"
    )


def aggregate_means(
    adata,
    genes: Iterable[str],
    group_by: str,
    *,
    layer: str | None = None,
    use_raw: bool | None = None,
    standard_scale: StandardScale | None = None,
    extra_by: Iterable[str] = (),
) -> pd.DataFrame:
    """Return long mean expression without computing detection fractions."""
    genes = ordered_unique(genes)
    by = list(dict.fromkeys([group_by, *extra_by]))
    means = group_means(
        adata,
        genes,
        group_by,
        layer=layer,
        use_raw=use_raw,
        extra_by=extra_by,
    )
    means = _standardize(means, standard_scale)
    n_groups = len(means)
    index_frame = means.index.to_frame(index=False)
    long = pd.DataFrame(
        {
            **{name: np.tile(index_frame[name].to_numpy(), len(genes)) for name in by},
            "feature": np.repeat(genes, n_groups),
            "mean_expression": means.to_numpy().T.reshape(-1),
        }
    )
    long["feature"] = pd.Categorical(long["feature"], categories=genes, ordered=True)
    return long


def aggregate_expression(
    adata,
    genes: Iterable[str],
    group_by: str,
    *,
    layer: str | None = None,
    use_raw: bool | None = None,
    expression_cutoff: float = 0.0,
    standard_scale: StandardScale | None = None,
    extra_by: Iterable[str] = (),
) -> pd.DataFrame:
    """Return a long DataFrame ``[group_by, *extra_by, feature, mean_expression, fraction]``.

    ``feature`` is an ordered categorical following the order of ``genes`` so
    downstream plots keep the requested gene order. ``extra_by`` adds grouping
    columns (e.g. a ``split_by`` facet) so means/fractions are computed per
    ``group_by`` × ``extra_by`` combination. ``standard_scale`` may be
    ``'var'`` / ``'group'`` (scanpy conventions) or ``'zscore'`` (an ggann
    extension); ``None`` leaves the raw group means untouched.
    """
    genes = ordered_unique(genes)
    by = list(dict.fromkeys([group_by, *extra_by]))  # dedupe, keep order
    kind, lyr = resolve_source(adata, layer, use_raw)
    mean_df, frac_df = _grouped_expression(
        adata,
        genes,
        by,
        kind=kind,
        layer=lyr,
        expression_cutoff=expression_cutoff,
    )
    assert frac_df is not None
    mean_df = _standardize(mean_df, standard_scale)

    n_groups = len(mean_df)
    index_frame = mean_df.index.to_frame(index=False)
    long = pd.DataFrame(
        {
            **{name: np.tile(index_frame[name].to_numpy(), len(genes)) for name in by},
            "feature": np.repeat(genes, n_groups),
            "mean_expression": mean_df.to_numpy().T.reshape(-1),
            "fraction": frac_df.to_numpy().T.reshape(-1),
        }
    )
    long["feature"] = pd.Categorical(long["feature"], categories=genes, ordered=True)
    return long
