"""Grouped expression aggregation, expressed through annplyr ``summarize``.

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

import annplyr as ap
import numpy as np
import pandas as pd

from ._expression import densify_frame as _densify
from ._expression import ordered_unique, project_expression, resolve_source

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


def _result_names(stem: str, count: int, reserved: Iterable[str]) -> list[str]:
    """Return deterministic temporary columns that cannot shadow group keys."""
    taken = set(reserved)
    names: list[str] = []
    index = 0
    while len(names) < count:
        candidate = f"__ggann_{stem}_{index}"
        index += 1
        if candidate in taken:
            continue
        names.append(candidate)
        taken.add(candidate)
    return names


def tidy_expression(adata, genes, group_by, *, layer=None, use_raw=None, extra_obs=()):
    """Long per-cell expression ``[obs_name, feature, value, group_by, *extra_obs]``.

    Shared by the violin / stacked-violin / tracksplot paths so the source
    picking, densification and feature-ordering live in one place. ``extra_obs``
    carries additional obs columns (e.g. a ``split_by`` facet) into the frame.
    """
    genes = ordered_unique(genes)
    obs = list(dict.fromkeys([group_by, *extra_obs]))  # dedupe, keep order
    kind, lyr = resolve_source(adata, layer, use_raw)
    projected, genes = project_expression(adata, genes, kind=kind, layer=lyr, obs=obs)
    tidy = projected.ap.to_tidy(obs=obs, x=genes)
    tidy = _densify(tidy)
    original_names = pd.Series(
        adata.obs_names.to_numpy(), index=projected.obs_names, dtype=object
    )
    tidy["obs_name"] = tidy["obs_name"].map(original_names).to_numpy()
    tidy["feature"] = pd.Categorical(tidy["feature"], categories=genes, ordered=True)
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
    projected, genes = project_expression(adata, genes, kind=kind, layer=lyr, obs=by)
    result_names = _result_names("mean", len(genes), by)
    mean_expr = {
        result: ap.mean(ap.col(gene))
        for result, gene in zip(result_names, genes, strict=True)
    }
    mean_df = projected.ap.summarize(x=mean_expr, by=by)
    mean_df = _densify(mean_df).set_index(by)[result_names].astype(float)
    mean_df.columns = genes
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
    long = means.reset_index().melt(
        id_vars=by, var_name="feature", value_name="mean_expression"
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
    projected, genes = project_expression(adata, genes, kind=kind, layer=lyr, obs=by)

    mean_columns = _result_names("mean", len(genes), by)
    fraction_columns = _result_names("fraction", len(genes), by)
    mean_summary = _densify(
        projected.ap.summarize(
            x={
                column: ap.mean(ap.col(gene))
                for column, gene in zip(mean_columns, genes, strict=True)
            },
            by=by,
        )
    )

    # annplyr treats ``col > cutoff`` as a complex group-by expression. Build a
    # bounded boolean projection instead, then use its native mean reduction.
    # For a negative cutoff, implicit sparse zeros are true, so only the already
    # gene-projected n_obs x n_requested matrix must become dense.
    from anndata import AnnData
    from scipy import sparse

    expression = projected.X
    if sparse.issparse(expression) and expression_cutoff >= 0:
        detected = (expression > expression_cutoff).astype(np.float32)
    else:
        values = (
            expression.toarray()
            if sparse.issparse(expression)
            else np.asarray(expression)
        )
        detected = np.greater(values, expression_cutoff).astype(np.float32)
    detected_adata = AnnData(X=detected, obs=projected.obs, var=projected.var)
    fraction_summary = _densify(
        detected_adata.ap.summarize(
            x={
                column: ap.mean(ap.col(gene))
                for column, gene in zip(fraction_columns, genes, strict=True)
            },
            by=by,
        )
    )

    mean_df = mean_summary.set_index(by)[mean_columns].astype(float)
    mean_df.columns = genes
    frac_df = fraction_summary.set_index(by)[fraction_columns].astype(float)
    frac_df = frac_df.reindex(mean_df.index)
    frac_df.columns = genes
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
