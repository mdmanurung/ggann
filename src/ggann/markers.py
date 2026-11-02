"""Marker-summary plots: tracksplot and the gene-group dot/matrix plots.

``plot_tracksplot`` stacks genes as rows (``facet_grid`` over the feature),
mirroring ``sc.pl.tracksplot``. The ``*_grouped`` helpers split the gene axis
into one labelled facet panel per group, where scanpy draws
``var_group_labels`` brackets. Expression is pulled long
through annplyr's ``to_tidy``. (The violin family lives in ``distributions``.)
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
    scale_x_discrete,
    theme,
)

from ._aggregate import tidy_expression
from ._compat import accepts_groupby
from ._expression import ordered_unique
from ._grouping import _drop_missing_groups, _group_categories, _order_groups
from ._palette import scale_fill_obs
from .plots import (
    _cell_rank,
    plot_dotplot,
    plot_matrixplot,
)
from .publication import _family_theme

__all__ = [
    "plot_tracksplot",
    "plot_dotplot_grouped",
    "plot_matrixplot_grouped",
]


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


def _add_gene_group_facet(plot, gene_groups: dict):
    """Split the gene axis into one labelled facet panel per group (``gene_group``).

    ``gene_groups`` maps a label -> list of genes; the panels appear in dict order.
    A gene listed under several labels is drawn once in each panel, labelled by its
    group, with values from the plot's single per-gene computation.
    """
    genes, occ = _flatten_gene_groups(gene_groups)
    if len(genes) == len(occ):  # no repeats: label the existing rows in place
        label = {g: lab for lab, g in occ}
        data = plot.data.copy()
        data["gene_group"] = pd.Categorical(
            data["feature"].astype(str).map(label),
            categories=list(dict.fromkeys(gene_groups)),
            ordered=True,
        )
        plot.data = data
    else:  # repeats: one keyed copy of the gene's rows per occurrence
        plot.data = _expand_occurrences(plot.data, occ, label_col="gene_group")
        plot = plot + scale_x_discrete(labels=lambda ks: [str(k).split(_SEP, 1)[-1] for k in ks])
    return plot + facet_grid(". ~ gene_group", scales="free_x", space="free_x")


@accepts_groupby
def plot_dotplot_grouped(adata, gene_groups: dict, group_by: str, **kwargs):
    """Draw a dotplot with one facet panel per gene group (not Scanpy brackets).

    Parameters
    ----------
    adata
        Annotated data matrix.
    gene_groups : mapping of str to sequence of str
        Labelled marker groups in display order.
    group_by : str
        Observation grouping column.
    **kwargs
        Passed to :func:`ggann.plot_dotplot`.

    Returns
    -------
    plotnine.ggplot
        Composable grouped marker dotplot.

    Raises
    ------
    KeyError
        If genes or grouping metadata are missing.
    ValueError
        If forwarded options are invalid.

    Notes
    -----
    A gene listed under several labels appears once per label (as in scanpy), with
    values from one computation per gene. Extraction and ownership match
    :func:`ggann.plot_dotplot`.

    Examples
    --------
    >>> p = plot_dotplot_grouped(adata, {"T": ["CD3D"]}, group_by="cell_type")
    """
    genes, _ = _flatten_gene_groups(gene_groups)
    return _add_gene_group_facet(plot_dotplot(adata, genes, group_by, **kwargs), gene_groups)


@accepts_groupby
def plot_matrixplot_grouped(adata, gene_groups: dict, group_by: str, **kwargs):
    """Draw a matrixplot with one facet panel per gene group (not Scanpy brackets).

    Parameters
    ----------
    adata
        Annotated data matrix.
    gene_groups : mapping of str to sequence of str
        Labelled marker groups in display order.
    group_by : str
        Observation grouping column.
    **kwargs
        Passed to :func:`ggann.plot_matrixplot`.

    Returns
    -------
    plotnine.ggplot
        Composable grouped marker matrixplot.

    Raises
    ------
    KeyError
        If genes or grouping metadata are missing.
    ValueError
        If forwarded options are invalid.

    Notes
    -----
    A gene listed under several labels appears once per label (as in scanpy), with
    values from one computation per gene. Extraction and ownership match
    :func:`ggann.plot_matrixplot`.

    Examples
    --------
    >>> p = plot_matrixplot_grouped(adata, {"T": ["CD3D"]}, group_by="cell_type")
    """
    genes, _ = _flatten_gene_groups(gene_groups)
    return _add_gene_group_facet(plot_matrixplot(adata, genes, group_by, **kwargs), gene_groups)


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
