"""Per-group expression distribution plots.

The whole distribution family lives here: violin (+ stacked), box, sina/beeswarm,
ridgeline, and the bar/line expression summaries. Inspired by DOTools / scplotter's
``FeatureStatPlot`` family and scanpy's ``sc.pl.violin`` / ``stacked_violin``. All
extraction goes through annplyr (``tidy_expression`` / ``resolve_frame``) and every
helper returns a plain :class:`plotnine.ggplot`, so ``+ scale_* / + theme(...) /
+ facet_*`` still compose.
"""

from __future__ import annotations

from typing import Iterable, Sequence

import numpy as np
import pandas as pd
import plotnine_extra as pe
from plotnine import (
    aes,
    element_blank,
    element_text,
    facet_grid,
    geom_boxplot,
    geom_col,
    geom_errorbar,
    geom_jitter,
    geom_line,
    geom_point,
    geom_ribbon,
    geom_violin,
    ggplot,
    labs,
    scale_y_continuous,
    theme,
)

from ._aggregate import tidy_expression
from ._expression import ordered_unique
from ._grouping import _downsample_cells, _group_categories, _order_groups
from ._palette import scale_color_obs, scale_fill_obs
from ._resolve import plain_name, resolve_frame
from .plots import (
    _feature_facet,
    _is_numeric,
)
from .theme import theme_ggann

__all__ = [
    "plot_violin",
    "plot_stacked_violin",
    "plot_box",
    "plot_sina",
    "plot_ridge",
    "plot_expression_bar",
    "plot_expression_line",
]


def _summarise(values: pd.core.groupby.SeriesGroupBy, error: str, agg="mean") -> pd.DataFrame:
    """Central tendency (``agg``) and a symmetric error (se / sd / none) per group.

    ``agg`` is any pandas reduction name or callable (``"mean"``, ``"median"``,
    ``np.sum`` ...); the result column is still called ``mean`` for the callers.
    """
    summary = values.agg(mean=agg, sd="std", n="count").reset_index()
    summary["sd"] = summary["sd"].fillna(0.0)
    if error == "se":
        summary["err"] = summary["sd"] / np.sqrt(summary["n"].clip(lower=1))
    elif error == "sd":
        summary["err"] = summary["sd"]
    elif error == "none":
        summary["err"] = 0.0
    else:
        raise ValueError(f"error must be 'se', 'sd' or 'none', got {error!r}")
    summary["ymin"] = summary["mean"] - summary["err"]
    summary["ymax"] = summary["mean"] + summary["err"]
    return summary


def plot_box(
    adata,
    genes: Sequence[str],
    group_by: str,
    *,
    split_by: str | None = None,
    layer: str | None = None,
    use_raw: bool | None = None,
    ncol: int = 1,
    jitter: bool = True,
    jitter_size: float = 0.35,
    jitter_alpha: float = 0.25,
    stats: bool = False,
    categories_order: Sequence[str] | None = None,
    downsample: int | None = None,
    random_state: int | None = 0,
):
    """Per-group expression box plots, one facet per gene, with jittered cells overlaid.

    Set ``jitter=False`` for a plain box plot, ``split_by`` for a gene x split facet
    grid, or ``stats=True`` to overlay a group-comparison test via plotnine-extra's
    ``stat_compare_means``. Pass ``downsample=N`` to cap cells per group before the
    (jitter) draw for large data — see :func:`ggann.plot_violin`.

    Note: ``downsample`` subsets the cells the geoms see, so with ``stats=True``
    the p-value is computed on the subsample, and the boxplot's outliers reflect
    only the kept cells. ``random_state`` controls which cells are retained.

    Parameters
    ----------
    adata
        Annotated data matrix.
    genes : sequence of str
        Genes shown in facets.
    group_by : str
        Observation grouping column.
    split_by : str, optional
        Additional facet column.
    layer, use_raw : optional
        Mutually exclusive expression source.
    ncol : int
        Facet columns.
    jitter, stats : bool
        Add cell points or a group-comparison layer.
    jitter_size, jitter_alpha : float
        Point size and opacity.
    categories_order : sequence of str, optional
        Complete observed group order.
    downsample : int, optional
        Maximum cells per group.
    random_state : int, optional
        Downsampling seed.

    Returns
    -------
    plotnine.ggplot
        Composable faceted box plot.

    Raises
    ------
    KeyError
        If requested genes, metadata, or a layer are missing.
    ValueError
        If source, ordering, or downsampling is invalid.

    Notes
    -----
    Only requested genes are projected; the long table scales with retained cells.
    The input is not mutated.

    Examples
    --------
    >>> p = plot_box(adata, ["CD3D"], group_by="cell_type", jitter=False)
    """
    adata = _downsample_cells(adata, group_by, downsample, random_state=random_state)
    genes = ordered_unique(genes)
    extra = [split_by] if split_by else []
    tidy = tidy_expression(adata, genes, group_by, layer=layer, use_raw=use_raw, extra_obs=extra)
    if categories_order is None:
        categories_order = _group_categories(adata, group_by)
    tidy = _order_groups(tidy, group_by, categories_order)

    plot = ggplot(tidy, aes(group_by, "value", fill=group_by)) + geom_boxplot(
        width=0.7, outlier_alpha=0.0 if jitter else 1.0
    )
    if jitter:
        plot = plot + geom_jitter(
            width=0.2, height=0.0, size=jitter_size, alpha=jitter_alpha, stroke=0
        )
    plot = (
        plot
        + _feature_facet(split_by, ncol=ncol)
        + scale_fill_obs(adata, group_by)
        + labs(x="", y="expression", fill=group_by)
        + theme_ggann()
        + pe.rotate_x_text(45)
    )
    if stats:
        plot = plot + pe.stat_compare_means()
    return plot


def plot_sina(
    adata,
    genes: Sequence[str],
    group_by: str,
    *,
    split_by: str | None = None,
    layer: str | None = None,
    use_raw: bool | None = None,
    ncol: int = 1,
    size: float = 0.4,
    alpha: float = 0.5,
    violin: bool = True,
    bins: int = 50,
    categories_order: Sequence[str] | None = None,
    downsample: int | None = None,
    random_state: int | None = 0,
):
    """Sina / beeswarm of per-group expression, one facet per gene.

    A sina plot spreads each group's cells horizontally in proportion to the local
    density (via plotnine-extra's ``geom_sina``), so it shows every cell like a
    jitter but keeps the shape of a violin. ``violin=True`` draws a faint violin
    behind the points for context. ``downsample=N`` caps cells per group first --
    useful for large data, since a sina draws one mark per cell. ``random_state``
    controls which cells are retained.

    Parameters
    ----------
    adata
        Annotated data matrix.
    genes : sequence of str
        Genes shown in facets.
    group_by : str
        Observation grouping column.
    split_by : str, optional
        Additional facet column.
    layer, use_raw : optional
        Mutually exclusive expression source.
    ncol : int
        Facet columns.
    size, alpha : float
        Point size and opacity.
    violin : bool
        Draw a contextual violin.
    bins : int
        Density bins used by the sina geometry.
    categories_order : sequence of str, optional
        Complete observed group order.
    downsample : int, optional
        Maximum cells per group.
    random_state : int, optional
        Downsampling seed.

    Returns
    -------
    plotnine.ggplot
        Composable sina plot.

    Raises
    ------
    KeyError
        If requested genes, metadata, or a layer are missing.
    ValueError
        If source, ordering, or downsampling is invalid.

    Notes
    -----
    One point is prepared per retained cell and gene; ``adata`` is not mutated.

    Examples
    --------
    >>> p = plot_sina(adata, ["CD3D"], group_by="cell_type")
    """
    genes = ordered_unique(genes)
    adata = _downsample_cells(adata, group_by, downsample, random_state=random_state)
    extra = [split_by] if split_by else []
    tidy = tidy_expression(adata, genes, group_by, layer=layer, use_raw=use_raw, extra_obs=extra)
    if categories_order is None:
        categories_order = _group_categories(adata, group_by)
    tidy = _order_groups(tidy, group_by, categories_order)

    # Drive stat_sina by binwidth rather than bins: plotnine 0.15's stat_sina
    # mutates its shared ``bins`` param across panels (the second facet then gets
    # an array where a scalar is expected and crashes). The binwidth branch reads
    # an untouched param each panel, so it is multi-facet safe.
    #
    # Size the binwidth off the NARROWEST gene's range (÷ bins) so every facet
    # gets ~``bins`` bins -- a single global span sized off the widest gene would
    # collapse a narrow-range gene's panel to one or two bins. Floor it by the
    # widest gene's range so a pathologically narrow gene can't drive the widest
    # panel to an unbounded bin count (which would be slow to build).
    ranges = tidy.groupby("feature", observed=True)["value"].agg(lambda s: s.max() - s.min())
    ranges = ranges[ranges > 0]
    if len(ranges):
        binwidth = max(float(ranges.min()) / bins, float(ranges.max()) / 1000)
    else:
        binwidth = 1.0

    plot = ggplot(tidy, aes(group_by, "value", color=group_by))
    if violin:
        plot = plot + geom_violin(
            aes(fill=group_by),
            color="none",
            alpha=0.15,
            scale="width",
            show_legend=False,
        )
    plot = (
        plot
        + pe.geom_sina(size=size, alpha=alpha, stroke=0, binwidth=binwidth)
        + _feature_facet(split_by, ncol=ncol)
        + scale_color_obs(adata, group_by)
        + scale_fill_obs(adata, group_by)
        + labs(x="", y="expression", color=group_by)
        + theme_ggann()
        + pe.rotate_x_text(45)
    )
    return plot


def plot_expression_bar(
    adata,
    genes: Sequence[str],
    group_by: str,
    *,
    split_by: str | None = None,
    layer: str | None = None,
    use_raw: bool | None = None,
    ncol: int = 1,
    agg="mean",
    error: str = "se",
    categories_order: Sequence[str] | None = None,
):
    """Aggregated expression per group as bars with an error bar, one facet per gene.

    ``agg`` is the bar height (``"mean"`` default, or ``"median"`` / ``np.sum`` / any
    pandas reduction). ``error`` is the summary bar: ``"se"`` (default), ``"sd"`` or
    ``"none"``. ``split_by`` adds a gene x split facet grid. Bars start at zero and
    hide the distribution -- for a distribution-honest view use :func:`plot_box` or
    :func:`ggann.plot_violin`.

    Parameters
    ----------
    adata
        Annotated data matrix.
    genes : sequence of str
        Genes shown in facets.
    group_by : str
        Observation grouping column.
    split_by : str, optional
        Additional grouping/facet column.
    layer, use_raw : optional
        Mutually exclusive expression source.
    ncol : int
        Facet columns.
    agg
        Pandas aggregation name or callable.
    error : {"se", "sd", "none"}
        Error-bar summary.
    categories_order : sequence of str, optional
        Complete observed group order.

    Returns
    -------
    plotnine.ggplot
        Composable expression bar plot.

    Raises
    ------
    KeyError
        If requested genes, metadata, or a layer are missing.
    ValueError
        If source, ordering, aggregation, or error mode is invalid.

    Notes
    -----
    Requested genes are projected before aggregation; ``adata`` is unchanged.

    Examples
    --------
    >>> p = plot_expression_bar(adata, ["CD3D"], group_by="cell_type")
    """
    genes = ordered_unique(genes)
    extra = [split_by] if split_by else []
    tidy = tidy_expression(adata, genes, group_by, layer=layer, use_raw=use_raw, extra_obs=extra)
    if categories_order is None:
        categories_order = _group_categories(adata, group_by)
    tidy = _order_groups(tidy, group_by, categories_order)

    by = [group_by, *extra, "feature"]
    summary = _summarise(tidy.groupby(by, observed=True)["value"], error, agg=agg)
    summary["feature"] = pd.Categorical(summary["feature"], categories=genes, ordered=True)

    ylab = f"{agg} expression" if isinstance(agg, str) else "expression"
    plot = (
        ggplot(summary, aes(group_by, "mean", fill=group_by))
        + geom_col(width=0.7)
        + _feature_facet(split_by, ncol=ncol)
        + scale_fill_obs(adata, group_by)
        + labs(x="", y=ylab, fill=group_by)
        + theme_ggann()
        + pe.rotate_x_text(45)
    )
    if error != "none":
        plot = plot + geom_errorbar(aes(ymin="ymin", ymax="ymax"), width=0.3)
    return plot


def plot_expression_line(
    adata,
    genes: Sequence[str],
    x: str,
    *,
    group_by: str | None = None,
    layer: str | None = None,
    use_raw: bool | None = None,
    ncol: int = 1,
    agg="mean",
    error: str = "se",
    categories_order: Sequence[str] | None = None,
):
    """Aggregated expression trend across an ordered variable ``x``, one facet per gene.

    ``x`` is an obs column (an ordered condition, timepoint or a numeric covariate
    such as pseudotime). With ``group_by`` given, one line is drawn per group and
    coloured by it (e.g. an expression trajectory per cell type across timepoints).
    ``agg`` sets the summarised value (``"mean"`` default, or ``"median"`` etc.);
    ``error`` adds a summary error bar per point (``"se"`` / ``"sd"`` / ``"none"``).

    Parameters
    ----------
    adata
        Annotated data matrix.
    genes : sequence of str
        Genes shown in facets.
    x : str
        Ordered or numeric observation column.
    group_by : str, optional
        Observation column defining lines.
    layer, use_raw : optional
        Mutually exclusive expression source.
    ncol : int
        Facet columns.
    agg
        Pandas aggregation name or callable.
    error : {"se", "sd", "none"}
        Error-bar summary.
    categories_order : sequence of str, optional
        Complete observed group order.

    Returns
    -------
    plotnine.ggplot
        Composable expression trend plot.

    Raises
    ------
    KeyError
        If genes or observation columns are missing.
    ValueError
        If source, ordering, aggregation, or error mode is invalid.

    Notes
    -----
    Requested genes and metadata are projected before aggregation; input is unchanged.

    Examples
    --------
    >>> p = plot_expression_line(adata, ["CD3D"], x="time", group_by="cell_type")
    """
    genes = ordered_unique(genes)
    xname = plain_name(adata, x)
    gname = plain_name(adata, group_by) if group_by is not None else None

    gene_names = [plain_name(adata, g) for g in genes]
    id_vars = [xname] + ([gname] if gname is not None else [])
    clash = set(gene_names) & set(id_vars)
    if clash:
        raise ValueError(
            f"feature name(s) {sorted(clash)} collide with x/group_by; melt would drop them. "
            f"Rename, or disambiguate with gene()/obs()."
        )

    cols = [x] + ([group_by] if group_by is not None else []) + list(genes)
    frame = resolve_frame(adata, cols, layer=layer, use_raw=use_raw)  # already densified
    long = frame.melt(
        id_vars=id_vars, value_vars=gene_names, var_name="feature", value_name="value"
    )
    long["feature"] = pd.Categorical(long["feature"], categories=gene_names, ordered=True)

    # Order the x axis: numeric stays numeric; categorical keeps its obs order.
    if not _is_numeric(long[xname]):
        x_cats = _group_categories(adata, x)
        if x_cats is None:
            x_cats = list(pd.unique(long[xname]))
        long[xname] = pd.Categorical(long[xname], categories=x_cats, ordered=True)
    if gname is not None and categories_order is None:
        categories_order = _group_categories(adata, group_by)
    if gname is not None:
        long = _order_groups(long, gname, categories_order)

    summary = _summarise(
        long.groupby(id_vars + ["feature"], observed=True)["value"], error, agg=agg
    )
    summary["feature"] = pd.Categorical(summary["feature"], categories=gene_names, ordered=True)

    color = gname if gname is not None else None
    mapping = (
        aes(x=xname, y="mean", color=color, group=color)
        if color
        else aes(x=xname, y="mean", group=1)
    )
    ylab = f"{agg} expression" if isinstance(agg, str) else "expression"
    plot = (
        ggplot(summary, mapping)
        + geom_line()
        + geom_point(size=2)
        + pe.facet_wrap("~feature", ncol=ncol, scales="free_y")
        + labs(x=xname, y=ylab, color=gname)
        + theme_ggann()
    )
    if color is not None:
        plot = plot + scale_color_obs(adata, gname)
    if error != "none":
        plot = plot + geom_errorbar(aes(ymin="ymin", ymax="ymax"), width=0.15)
    return plot


def plot_violin(
    adata,
    genes: Sequence[str],
    group_by: str,
    *,
    split_by: str | None = None,
    layer: str | None = None,
    use_raw: bool | None = None,
    ncol: int = 1,
    scale: str = "width",
    add_box: bool = True,
    add_points: bool = False,
    stats: bool = False,
    downsample: int | None = None,
    random_state: int | None = 0,
    categories_order: Iterable[str] | None = None,
):
    """Per-group expression distributions, one facet per gene (stacked-violin style).

    ``add_box=True`` (default) nests a slim white boxplot inside each violin so the
    median and quartiles read off cleanly, the way scplotter's ``FeatureStatPlot``
    does. ``add_points=True`` overlays the individual cells as jitter (scplotter's
    ``add_point``). ``split_by`` adds a facet column (gene rows x split columns).
    Set ``stats=True`` to overlay a group-comparison test via plotnine-extra's
    ``stat_compare_means``. ``downsample`` caps cells per group before the (slow)
    violin KDE -- a big speed-up on large data for a visually identical plot.

    Note: ``downsample`` subsets the cells the geoms see, so any ``stats=True``
    p-value is then computed on the *subsample*, not the full data. Leave
    ``downsample`` unset when you need the reported test to reflect every cell.
    ``random_state`` controls which cells are retained.

    Parameters
    ----------
    adata
        Annotated data matrix.
    genes : sequence of str
        Genes shown in facets.
    group_by : str
        Observation grouping column.
    split_by : str, optional
        Additional facet column.
    layer, use_raw : optional
        Mutually exclusive expression source.
    ncol : int
        Facet columns.
    scale : str
        plotnine violin-width scaling mode.
    add_box, add_points, stats : bool
        Add inner boxes, cell points, or group tests.
    downsample : int, optional
        Maximum cells per group.
    random_state : int, optional
        Downsampling seed.
    categories_order : iterable of str, optional
        Complete observed group order.

    Returns
    -------
    plotnine.ggplot
        Composable faceted violin plot.

    Raises
    ------
    KeyError
        If genes, grouping metadata, or a layer are missing.
    ValueError
        If source, ordering, scaling, or downsampling is invalid.

    Notes
    -----
    KDE and optional statistics use retained cells. ``adata`` is not mutated.

    Examples
    --------
    >>> p = plot_violin(adata, ["CD3D"], group_by="cell_type")
    """
    adata = _downsample_cells(adata, group_by, downsample, random_state=random_state)
    genes = ordered_unique(genes)
    extra = [split_by] if split_by else []
    tidy = tidy_expression(adata, genes, group_by, layer=layer, use_raw=use_raw, extra_obs=extra)
    if categories_order is None:
        categories_order = _group_categories(adata, group_by)
    tidy = _order_groups(tidy, group_by, categories_order)
    plot = ggplot(tidy, aes(group_by, "value", fill=group_by)) + geom_violin(scale=scale)
    # box first, then points on top -- otherwise the white box occludes the jitter
    if add_box:
        plot = plot + geom_boxplot(width=0.12, fill="white", outlier_alpha=0.0, show_legend=False)
    if add_points:
        plot = plot + geom_jitter(width=0.2, height=0.0, size=0.3, alpha=0.25, stroke=0)
    plot = (
        plot
        + _feature_facet(split_by, ncol=ncol)
        + scale_fill_obs(adata, group_by)
        + labs(x="", y="expression", fill=group_by)
        + theme_ggann()
        + pe.rotate_x_text(45)
    )
    if stats:
        plot = plot + pe.stat_compare_means()
    return plot


def plot_stacked_violin(
    adata,
    genes: Sequence[str],
    group_by: str,
    *,
    layer: str | None = None,
    use_raw: bool | None = None,
    scale: str = "width",
    categories_order=None,
    downsample: int | None = None,
    random_state: int | None = 0,
):
    """Compact genes-as-rows violin grid across groups (``sc.pl.stacked_violin``).

    ``downsample`` caps cells per group before the KDE. ``random_state`` controls
    which cells are retained.

    Parameters
    ----------
    adata
        Annotated data matrix.
    genes : sequence of str
        Genes shown as rows.
    group_by : str
        Observation grouping column.
    layer, use_raw : optional
        Mutually exclusive expression source.
    scale : str
        plotnine violin-width scaling mode.
    categories_order : iterable of str, optional
        Complete observed group order.
    downsample : int, optional
        Maximum cells per group.
    random_state : int, optional
        Downsampling seed.

    Returns
    -------
    plotnine.ggplot
        Composable stacked violin plot.

    Raises
    ------
    KeyError
        If genes, grouping metadata, or a layer are missing.
    ValueError
        If source, ordering, scaling, or downsampling is invalid.

    Notes
    -----
    KDE uses retained cells; ``adata`` is unchanged.

    Examples
    --------
    >>> p = plot_stacked_violin(adata, ["CD3D", "NKG7"], group_by="cell_type")
    """
    adata = _downsample_cells(adata, group_by, downsample, random_state=random_state)
    genes = ordered_unique(genes)
    tidy = tidy_expression(adata, genes, group_by, layer=layer, use_raw=use_raw)
    if categories_order is None:
        categories_order = _group_categories(adata, group_by)
    tidy = _order_groups(tidy, group_by, categories_order)
    return (
        ggplot(tidy, aes(group_by, "value", fill=group_by))
        + geom_violin(scale=scale)
        + facet_grid("feature ~ .", scales="free_y")
        + scale_fill_obs(adata, group_by)
        + labs(x="", y="", fill=group_by)
        + theme_ggann()
        + theme(strip_text_y=element_text(angle=0))
        + pe.rotate_x_text(45)
    )


def _kde_curve(values: np.ndarray, grid: np.ndarray) -> np.ndarray:
    """KDE density on ``grid``; zeros for degenerate (constant / tiny) groups."""
    from scipy.stats import gaussian_kde

    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size < 2 or np.ptp(values) == 0:
        return np.zeros_like(grid)
    try:
        dens = gaussian_kde(values)(grid)
    except np.linalg.LinAlgError:  # singular covariance
        return np.zeros_like(grid)
    peak = dens.max()
    return dens / peak if peak > 0 else dens


def plot_ridge(
    adata,
    genes: Sequence[str],
    group_by: str,
    *,
    layer: str | None = None,
    use_raw: bool | None = None,
    ncol: int = 1,
    scale: float = 1.6,
    n_grid: int = 256,
    categories_order: Sequence[str] | None = None,
):
    """Ridgeline plot: one density ridge per group, stacked and overlapping, per gene.

    ``scale`` sets ridge height relative to row spacing; values above one overlap
    neighbouring ridges. Groups with fewer than two cells or zero
    variance draw a flat baseline. plotnine has no ridgeline geom, so this builds one
    from a per-group gaussian KDE offset vertically and drawn as ``geom_ribbon``.

    Parameters
    ----------
    adata
        Annotated data matrix.
    genes : sequence of str
        Genes shown in facets.
    group_by : str
        Observation grouping column.
    layer, use_raw : optional
        Mutually exclusive expression source.
    ncol : int
        Facet columns.
    scale : float
        Ridge height relative to row spacing.
    n_grid : int
        KDE evaluation-grid size.
    categories_order : sequence of str, optional
        Complete observed group order.

    Returns
    -------
    plotnine.ggplot
        Composable ridgeline plot.

    Raises
    ------
    KeyError
        If genes, grouping metadata, or a layer are missing.
    ValueError
        If source, ordering, or KDE configuration is invalid.

    Notes
    -----
    KDE work scales with genes, groups, cells, and ``n_grid``. Input is unchanged.

    Examples
    --------
    >>> p = plot_ridge(adata, ["CD3D"], group_by="cell_type")
    """
    genes = ordered_unique(genes)
    tidy = tidy_expression(adata, genes, group_by, layer=layer, use_raw=use_raw)
    if categories_order is None:
        categories_order = _group_categories(adata, group_by)
    tidy = _order_groups(tidy, group_by, categories_order)
    order = list(tidy[group_by].cat.categories)
    pos = {g: i for i, g in enumerate(order)}

    rows = []
    for feature, fdf in tidy.groupby("feature", observed=True):
        lo, hi = float(fdf["value"].min()), float(fdf["value"].max())
        if hi <= lo:
            hi = lo + 1.0
        grid = np.linspace(lo, hi, n_grid)
        for g, gdf in fdf.groupby(group_by, observed=True):
            dens = _kde_curve(gdf["value"].to_numpy(), grid) * scale
            base = float(pos[g])
            rows.append(
                pd.DataFrame(
                    {
                        "x": grid,
                        "ymin": base,
                        "ymax": base + dens,
                        group_by: g,
                        "feature": str(feature),
                    }
                )
            )
    long = pd.concat(rows, ignore_index=True)
    long[group_by] = pd.Categorical(long[group_by], categories=order, ordered=True)
    long["feature"] = pd.Categorical(
        long["feature"], categories=[str(x) for x in genes], ordered=True
    )

    # Draw top ridges first so lower ones layer IN FRONT and overlap cleanly
    # (a single ribbon layer paints higher rows on top, clipping the peaks below).
    # Each ridge is a filled area with a thin white line tracing only its top edge,
    # not a boxed outline -- so the silhouettes blend instead of looking cut.
    plot = ggplot(long, aes("x"))
    for g in reversed(order):
        gd = long[long[group_by] == g]
        plot = plot + geom_ribbon(aes(ymin="ymin", ymax="ymax", fill=group_by), data=gd, alpha=0.9)
        plot = plot + geom_line(aes(y="ymax"), data=gd, color="white", size=0.4)
    return (
        plot
        + pe.facet_wrap("~feature", ncol=ncol, scales="free_x")
        + scale_fill_obs(adata, group_by)
        + scale_y_continuous(breaks=list(pos.values()), labels=order)
        + labs(x="expression", y="", fill=group_by)
        + theme_ggann()
        + theme(panel_grid=element_blank())
    )
