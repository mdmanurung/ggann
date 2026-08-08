"""Cell-type composition / proportion bar plots.

Counts cells per group (optionally within a sample/condition) through annplyr's
``count`` verb, then renders stacked / filled / dodged bars -- the standard
"what fraction of each sample is each cell type" figure that scanpy's basics
don't cover.
"""

from __future__ import annotations

import pandas as pd
import plotnine_extra as pe
from plotnine import aes, geom_area, geom_col, geom_line, geom_point, ggplot, labs

from ._grouping import _group_categories, _order_groups, _resolve_group_order
from ._palette import scale_color_obs, scale_fill_obs
from .theme import theme_ggann

__all__ = ["plot_proportions"]


def plot_proportions(
    adata,
    group_by: str,
    split_by: str | None = None,
    *,
    kind: str = "bar",
    normalize: bool = True,
    position: str = "stack",
    categories_order=None,
    split_order=None,
):
    """Cell-type composition across an optional ``split_by`` axis.

    Parameters
    ----------
    group_by
        The categorical whose composition is shown (the bar segments / fill).
    split_by
        Optional obs column for the x-axis (e.g. sample or condition); when
        omitted, one bar per ``group_by`` category.
    kind
        ``"bar"`` (default), ``"area"`` (stacked area across ``split_by``) or
        ``"trend"`` (a line per group across ``split_by``). ``"area"``/``"trend"``
        require ``split_by``. (Pie/donut/radar are not available -- plotnine has
        no polar coordinate system.)
    normalize
        Plot within-``split_by`` proportions (summing to 1) instead of raw counts.
    position
        Bar stacking: ``"stack"`` (default), ``"fill"`` (100 %), or ``"dodge"``.
    categories_order, split_order
        Explicit orders for the group and split columns. ``None`` uses the order
        stored on categorical ``obs`` columns, then observed values.
    """
    if kind not in ("bar", "area", "trend"):
        raise ValueError(f"kind must be 'bar', 'area' or 'trend', got {kind!r}")
    if kind in ("area", "trend") and split_by is None:
        raise ValueError(f"kind={kind!r} needs split_by (the x-axis to trend across).")
    grouping_columns = [group_by, *([split_by] if split_by is not None else [])]
    missing = [column for column in grouping_columns if column not in adata.obs]
    if missing:
        raise KeyError(f"Observation column(s) not found: {missing}.")
    by = [split_by, group_by] if split_by else [group_by]
    counts = adata.ap.count(by=by)
    if categories_order is None:
        categories_order = _group_categories(adata, group_by)

    if split_by:
        # count() omits (split x group) combinations with zero cells; fill them
        # with 0 so a stacked area/bar tiles cleanly instead of leaving gaps.
        if split_order is None:
            split_order = _group_categories(adata, split_by)
        gcats = _resolve_group_order(counts, group_by, categories_order)
        scats = _resolve_group_order(counts, split_by, split_order)
        observed_groups = set(counts[group_by].dropna())
        observed_splits = set(counts[split_by].dropna())
        gcats = [category for category in gcats if category in observed_groups]
        scats = [category for category in scats if category in observed_splits]
        categories_order, split_order = gcats, scats
        full = pd.MultiIndex.from_product([scats, gcats], names=[split_by, group_by])
        counts = (
            counts.set_index([split_by, group_by])["n"]
            .reindex(full, fill_value=0)
            .reset_index()
        )

    if normalize:
        if split_by:
            denom = counts.groupby(split_by, observed=True)["n"].transform("sum")
        else:
            denom = counts["n"].sum()
        counts["value"] = counts["n"] / denom
        ylab = "proportion of cells"
    else:
        counts["value"] = counts["n"]
        # geom_col(position="fill") re-normalizes to 0..1 regardless of the counts
        ylab = "proportion of cells" if position == "fill" else "number of cells"

    counts = _order_groups(counts, group_by, categories_order)
    if split_by:
        counts = _order_groups(counts, split_by, split_order)

    xvar = split_by if split_by else group_by

    if kind == "trend":
        return (
            ggplot(counts, aes(xvar, "value", color=group_by, group=group_by))
            + geom_line()
            + geom_point(size=2)
            + scale_color_obs(adata, group_by)
            + labs(x="", y=ylab, color=group_by)
            + theme_ggann()
            + pe.rotate_x_text(45)
        )
    if kind == "area":
        # no segment outline: on a many-group stacked area the thin white borders
        # read as slivers/gaps; the qualitative fill colours separate the bands.
        return (
            ggplot(counts, aes(xvar, "value", fill=group_by, group=group_by))
            + geom_area(position=position)
            + scale_fill_obs(adata, group_by)
            + labs(x="", y=ylab, fill=group_by)
            + theme_ggann()
            + pe.rotate_x_text(45)
        )
    # thin white borders separate the stacked segments cleanly (scplotter-style)
    col_kw = {"color": "white", "size": 0.15} if position != "dodge" else {}
    return (
        ggplot(counts, aes(xvar, "value", fill=group_by))
        + geom_col(position=position, width=0.9, **col_kw)
        + scale_fill_obs(adata, group_by)
        + labs(x="", y=ylab, fill=group_by)
        + theme_ggann()
        + pe.rotate_x_text(45)
    )
