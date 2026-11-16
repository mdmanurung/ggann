"""Marker-summary plots: tracksplot.

``plot_tracksplot`` stacks genes as rows (``facet_grid`` over the feature),
mirroring ``sc.pl.tracksplot``. Grouped marker dot and matrix plots take a
``{label: [genes]}`` dict in :func:`ggann.plot_dotplot` and
:func:`ggann.plot_matrixplot`. (The violin family lives in ``distributions``.)
"""

from __future__ import annotations

from typing import Sequence

import pandas as pd
from plotnine import (
    aes,
    element_blank,
    element_text,
    facet_grid,
    geom_rect,
    ggplot,
    labs,
    scale_x_continuous,
    theme,
)

from ._aggregate import tidy_expression
from ._compat import accepts_groupby
from ._expression import ordered_unique
from ._grouping import _drop_missing_groups, _group_categories, _order_groups
from ._palette import scale_fill_obs
from .plots import _cell_rank
from .publication import _family_theme

__all__ = ["plot_tracksplot"]


_SEP = "\x1f"  # joins a marker-group label and a gene into one panel-local key


def _occurrences(gene_groups: dict) -> list[tuple]:
    """``(label, gene)`` pairs in display order; a gene may recur under other labels."""
    if not isinstance(gene_groups, dict):
        raise TypeError(
            f"gene_groups must be a dict of {{label: [genes]}}, got {type(gene_groups).__name__}."
        )
    return [(label, g) for label, members in gene_groups.items() for g in ordered_unique(members)]


def _flatten_gene_groups(gene_groups: dict) -> tuple[list, list]:
    """Unique genes (first-seen order, one computation each) and the occurrence list."""
    occ = _occurrences(gene_groups)
    return ordered_unique(g for _, g in occ), occ


def _expand_occurrences(data: pd.DataFrame, occ: list, *, label_col: str) -> pd.DataFrame:
    """Repeat each gene's rows once per occurrence; ``feature`` becomes a per-occurrence key.

    Values are copied from the single per-gene computation, never recomputed. The
    plain gene name stays in the ``gene`` column.
    """
    keys = [f"{label}{_SEP}{g}" for label, g in occ]
    parts = []
    for (label, g), key in zip(occ, keys):
        part = data[data["feature"].astype(str) == str(g)].copy()
        part["gene"] = str(g)
        part[label_col] = label
        part["feature"] = key
        parts.append(part)
    out = pd.concat(parts, ignore_index=True)
    out["feature"] = pd.Categorical(out["feature"], categories=keys, ordered=True)
    out[label_col] = pd.Categorical(
        out[label_col], categories=list(dict.fromkeys(lab for lab, _ in occ)), ordered=True
    )
    return out


@accepts_groupby
def plot_tracksplot(
    adata,
    genes: Sequence[str] | dict,
    group_by: str,
    *,
    layer: str | None = None,
    use_raw: bool | None = None,
    categories_order=None,
):
    """Per-gene expression tracks across cells ordered by group.

    ``genes`` may be a ``{label: [genes]}`` dict of marker groups. Each gene is
    then one track per label it is listed under, titled ``gene [label]``, with
    values from one computation per gene. Group names label the x axis at each
    block's centre. ``categories_order`` overrides the order stored on a
    categorical ``obs`` column. ``None`` uses the stored order.

    Parameters
    ----------
    adata
        Annotated data matrix.
    genes : sequence of str or dict
        Genes shown as tracks, or marker groups ``{label: [genes]}``.
    group_by : str
        Observation column ordering cell blocks.
    layer, use_raw : optional
        Mutually exclusive expression source.
    categories_order : iterable of str, optional
        Complete observed group order.

    Returns
    -------
    plotnine.ggplot
        Composable per-cell tracks plot.

    Raises
    ------
    KeyError
        If genes, grouping metadata, or a layer are missing.
    ValueError
        If source selection or category ordering is invalid.

    Notes
    -----
    The long table has one row per cell and gene; no implicit downsampling occurs.
    Requested genes are projected and ``adata`` is unchanged. Tracks are drawn as
    one rectangle per cell and gene from zero to the value (negative values hang
    below zero).

    Examples
    --------
    >>> p = plot_tracksplot(adata, ["CD3D", "NKG7"], group_by="cell_type")
    """
    adata = _drop_missing_groups(adata, group_by)
    occ = None
    if isinstance(genes, dict):
        genes, occ = _flatten_gene_groups(genes)
    genes = ordered_unique(genes)
    tidy = tidy_expression(adata, genes, group_by, layer=layer, use_raw=use_raw)
    if categories_order is None:
        categories_order = _group_categories(adata, group_by)
    tidy = _order_groups(tidy, group_by, categories_order)

    # order cells by group so each track reads left-to-right by group
    tidy = _cell_rank(tidy, group_by)
    if occ is not None:
        tidy = _expand_occurrences(tidy, occ, label_col="gene_group")
        tidy["feature"] = pd.Categorical(
            [
                f"{g} [{k.split(_SEP, 1)[0]}]"
                for g, k in zip(tidy["gene"], tidy["feature"].astype(str))
            ],
            categories=[f"{g} [{label}]" for label, g in occ],
            ordered=True,
        )

    # one filled rectangle per cell and gene (geom_col is ~20x slower here)
    tidy["xmin"] = tidy["cell_rank"] - 0.5
    tidy["xmax"] = tidy["cell_rank"] + 0.5
    tidy["ymin"] = 0.0
    # group names at each block's centre
    centre = tidy.groupby(group_by, observed=True)["cell_rank"].agg(["min", "max"]).mean(axis=1)

    return (
        ggplot(tidy, aes(xmin="xmin", xmax="xmax", ymin="ymin", ymax="value", fill=group_by))
        + geom_rect()
        + facet_grid("feature ~ .", scales="free_y")
        + scale_fill_obs(adata, group_by)
        + scale_x_continuous(breaks=centre.to_numpy(), labels=[str(g) for g in centre.index])
        + labs(x="cells (ordered by group)", y="", fill=group_by)
        + _family_theme("matrix")
        + theme(
            axis_text_x=element_text(angle=90),
            axis_ticks_major_x=element_blank(),
            strip_text_y=element_text(angle=0),
        )
    )
