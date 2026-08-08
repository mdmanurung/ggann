"""Marker-summary plots: tracksplot and the gene-group-bracketed dot/matrix plots.

``plot_tracksplot`` stacks genes as rows (``facet_grid`` over the feature),
mirroring ``sc.pl.tracksplot``; the ``*_grouped`` helpers bracket the gene axis
into labelled groups (scanpy ``var_group_labels``). Expression is pulled long
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
    geom_col,
    ggplot,
    labs,
    theme,
)

from ._aggregate import tidy_expression
from ._expression import ordered_unique
from ._grouping import _group_categories, _order_groups
from ._palette import scale_fill_obs
from .plots import (
    _cell_rank,
    plot_dotplot,
    plot_matrixplot,
)
from .theme import theme_ggann

__all__ = [
    "plot_tracksplot",
    "plot_dotplot_grouped",
    "plot_matrixplot_grouped",
]


def _flatten_gene_groups(gene_groups: dict) -> tuple[list, dict]:
    """Flatten a ``{label: [genes]}`` mapping, rejecting a gene listed in two groups."""
    feat_to_group: dict = {}
    genes: list = []
    for label, members in gene_groups.items():
        for g in ordered_unique(members):
            if g in feat_to_group:
                raise ValueError(
                    f"gene {g!r} appears in more than one gene group "
                    f"({feat_to_group[g]!r} and {label!r}); each gene must belong to one group."
                )
            feat_to_group[g] = label
            genes.append(g)
    return genes, feat_to_group


def _add_gene_group_facet(plot, gene_groups: dict):
    """Bracket the gene axis into labelled groups via a ``gene_group`` facet column.

    ``gene_groups`` maps a label -> list of genes; the panels appear in dict order.
    Reuses the plot's own aggregated ``.data`` (which carries a ``feature`` column).
    """
    _, feat_to_group = _flatten_gene_groups(gene_groups)
    data = plot.data.copy()
    data["gene_group"] = pd.Categorical(
        data["feature"].astype(str).map(feat_to_group),
        categories=list(gene_groups),
        ordered=True,
    )
    plot.data = data
    return plot + facet_grid(". ~ gene_group", scales="free_x", space="free_x")


def plot_dotplot_grouped(adata, gene_groups: dict, group_by: str, **kwargs):
    """Dotplot with genes bracketed into labelled groups (scanpy ``var_group_labels``)."""
    genes, _ = _flatten_gene_groups(gene_groups)
    return _add_gene_group_facet(
        plot_dotplot(adata, genes, group_by, **kwargs), gene_groups
    )


def plot_matrixplot_grouped(adata, gene_groups: dict, group_by: str, **kwargs):
    """Matrixplot with genes bracketed into labelled groups (scanpy ``var_group_labels``)."""
    genes, _ = _flatten_gene_groups(gene_groups)
    return _add_gene_group_facet(
        plot_matrixplot(adata, genes, group_by, **kwargs), gene_groups
    )


def plot_tracksplot(
    adata,
    genes: Sequence[str],
    group_by: str,
    *,
    layer: str | None = None,
    use_raw: bool | None = None,
    categories_order=None,
):
    """Per-gene expression tracks across cells ordered by group.

    ``categories_order`` overrides the order stored on a categorical ``obs``
    column. ``None`` uses the stored order.
    """
    genes = ordered_unique(genes)
    tidy = tidy_expression(adata, genes, group_by, layer=layer, use_raw=use_raw)
    if categories_order is None:
        categories_order = _group_categories(adata, group_by)
    tidy = _order_groups(tidy, group_by, categories_order)

    # order cells by group so each track reads left-to-right by group
    tidy = _cell_rank(tidy, group_by)

    return (
        ggplot(tidy, aes("cell_rank", "value", fill=group_by))
        + geom_col(width=1.0)
        + facet_grid("feature ~ .", scales="free_y")
        + scale_fill_obs(adata, group_by)
        + labs(x="cells (ordered by group)", y="", fill=group_by)
        + theme_ggann()
        + theme(
            axis_text_x=element_blank(),
            axis_ticks_major_x=element_blank(),
            strip_text_y=element_text(angle=0),
        )
    )
