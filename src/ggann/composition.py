"""Cell-type composition / proportion plots.

Counts cells per group (optionally within a sample/condition) through annplyr's
``count`` verb, then renders stacked / filled / dodged bars, areas, trends,
pies, rings or radars -- the standard "what fraction of each sample is each
cell type" figure that scanpy's basics don't cover.
"""

from __future__ import annotations

import pandas as pd
import plotnine_extra as pe
from annplyr import data as _sources
from plotnine import (
    aes,
    coord_radial,
    element_blank,
    facet_wrap,
    geom_area,
    geom_col,
    geom_line,
    geom_point,
    ggplot,
    labs,
    theme,
)

from ._compat import accepts_scanpy_names
from ._grouping import _drop_missing_groups, _group_categories, _order_groups, _resolve_group_order
from ._palette import scale_color_obs, scale_fill_obs
from .publication import _family_theme
from .shapes.pie import geom_pie
from .shapes.radar import _radar_plot

__all__ = ["plot_proportions"]


@accepts_scanpy_names
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
    radar_norm: str | None = None,
):
    """Cell-type composition across an optional ``split_by`` axis.

    Parameters
    ----------
    adata
        Annotated data matrix.
    group_by
        The categorical whose composition is shown (the bar segments / fill).
    split_by
        Optional obs column for the x-axis (e.g. sample or condition); when
        omitted, one bar per ``group_by`` category.
    kind
        ``"bar"`` (default), ``"area"`` (stacked area across ``split_by``) or
        ``"trend"`` (a line per group across ``split_by``). ``"area"``/``"trend"``
        require ``split_by``. ``"pie"`` and ``"ring"`` draw one pie (or ring)
        of the group shares, one per ``split_by`` level as facets; their angles
        are always shares of each pie's total, so ``normalize`` does not change
        them. ``"radar"`` draws one spoke per ``group_by`` level and one polygon
        per ``split_by`` level (a single polygon without ``split_by``); it needs
        at least 3 ``group_by`` levels and a shared radius axis from 0.
    normalize
        Plot within-``split_by`` proportions (summing to 1) instead of raw counts.
    position
        Bar stacking: ``"stack"`` (default), ``"fill"`` (100 %), or ``"dodge"``.
        Ignored by ``"trend"``, ``"pie"``, ``"ring"`` and ``"radar"``.
    categories_order, split_order
        Explicit orders for the group and split columns. ``None`` uses the order
        stored on categorical ``obs`` columns, then observed values.
    radar_norm
        Only for ``kind="radar"`` with ``split_by`` and ``normalize=True``, where
        it is required. ``"within_split"`` draws each cell type's share of the
        cells in that ``split_by`` level, as the bars do. ``"across_splits"``
        draws the share of each cell type's cells that falls in that
        ``split_by`` level, so each spoke sums to 1 over the polygons. The two
        give different polygons, so ggann does not pick one for you.

    Returns
    -------
    plotnine.ggplot
        Composable composition plot.

    Raises
    ------
    KeyError
        If ``group_by`` or ``split_by`` is absent.
    ValueError
        If ``kind`` or ordering is invalid, area/trend lacks ``split_by``,
        ``split_by`` is the same column as ``group_by``, a radar has fewer than
        3 ``group_by`` levels, or ``radar_norm`` is missing or invalid for a
        normalised radar with ``split_by``.

    Notes
    -----
    Only observation metadata is counted; no expression matrix is materialized and
    ``adata`` is not mutated. Missing grouping values are omitted. Pie, ring and
    radar do not lump small groups into "Other", and the radar does not rescale
    spokes, so a rare cell type sits near the centre.

    Examples
    --------
    >>> p = plot_proportions(adata, group_by="cell_type", split_by="condition")
    """
    adata = _drop_missing_groups(adata, group_by, split_by)
    if split_by is not None and split_by == group_by:
        raise ValueError(
            f"split_by and group_by are both {group_by!r}; composition of a column "
            "within itself is trivially 1. Pass a different split_by."
        )
    kinds = ("bar", "area", "trend", "pie", "ring", "radar")
    if kind not in kinds:
        raise ValueError(f"kind must be one of {kinds}, got {kind!r}")
    radar_norms = ("within_split", "across_splits")
    if radar_norm is not None and radar_norm not in radar_norms:
        raise ValueError(f"radar_norm must be one of {radar_norms}, got {radar_norm!r}")
    if kind == "radar" and split_by and normalize and radar_norm is None:
        raise ValueError(
            "kind='radar' with split_by and normalize=True needs radar_norm: "
            "'within_split' (share of the cells in each split level, as the bars show) "
            "or 'across_splits' (share of each group's cells that falls in each split "
            "level). Or pass normalize=False for counts."
        )
    if kind in ("area", "trend") and split_by is None:
        raise ValueError(f"kind={kind!r} needs split_by (the x-axis to trend across).")
    grouping_columns = [group_by, *([split_by] if split_by is not None else [])]
    missing = [column for column in grouping_columns if column not in adata.obs]
    if missing:
        raise KeyError(f"Observation column(s) not found: {missing}.")
    by = [split_by, group_by] if split_by else [group_by]
    counts = adata.ap.count(*(_sources.obs[name] for name in by))
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
            counts.set_index([split_by, group_by])["n"].reindex(full, fill_value=0).reset_index()
        )

    if normalize:
        if split_by:
            across = kind == "radar" and radar_norm == "across_splits"
            denom = counts.groupby(group_by if across else split_by, observed=True)["n"].transform(
                "sum"
            )
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

    if kind in ("pie", "ring"):
        plot = (
            ggplot(counts, aes(fill=group_by, weight="n"))
            + geom_pie(color="white", size=0.15)
            + coord_radial(theta="y", expand=False, inner_radius=0.5 if kind == "ring" else 0.0)
            + scale_fill_obs(adata, group_by)
            + labs(fill=group_by)
            + _family_theme("standard")
            + theme(
                axis_text=element_blank(),
                axis_ticks=element_blank(),
                axis_title=element_blank(),
                panel_grid=element_blank(),
                panel_border=element_blank(),
            )
        )
        return plot + facet_wrap(f"~{split_by}") if split_by else plot
    if kind == "radar":
        n_axes = int(counts[group_by].nunique())
        if n_axes < 3:
            raise ValueError(f"kind='radar' needs at least 3 {group_by!r} levels, got {n_axes}")
        plot = _radar_plot(
            counts,
            group_by,
            "value",
            split_by,
            n_axes,
            (0.0, float(counts["value"].max())),
            0.25,
            True,
        )
        if split_by:
            plot = plot + scale_color_obs(adata, split_by) + scale_fill_obs(adata, split_by)
        return plot

    if kind == "trend":
        return (
            ggplot(counts, aes(xvar, "value", color=group_by, group=group_by))
            + geom_line()
            + geom_point(size=2)
            + scale_color_obs(adata, group_by)
            + labs(x="", y=ylab, color=group_by)
            + _family_theme("standard")
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
            + _family_theme("standard")
            + pe.rotate_x_text(45)
        )
    # thin white borders separate the stacked segments cleanly (scplotter-style)
    col_kw = {"color": "white", "size": 0.15} if position != "dodge" else {}
    return (
        ggplot(counts, aes(xvar, "value", fill=group_by))
        + geom_col(position=position, width=0.9, **col_kw)
        + scale_fill_obs(adata, group_by)
        + labs(x="", y=ylab, fill=group_by)
        + _family_theme("standard")
        + pe.rotate_x_text(45)
    )
