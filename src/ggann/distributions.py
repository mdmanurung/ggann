"""Per-group expression distribution plots.

The whole distribution family lives here: violin (+ stacked), box, sina/beeswarm,
ridgeline, and the bar/line expression summaries. Inspired by DOTools / scplotter's
``FeatureStatPlot`` family and scanpy's ``sc.pl.violin`` / ``stacked_violin``. All
extraction goes through annplyr (``tidy_expression`` / ``resolve_frame``) and every
helper returns a plain :class:`plotnine.ggplot`, so ``+ scale_* / + theme(...) /
+ facet_*`` still compose.
"""

from __future__ import annotations

import math
import warnings
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
    guide_legend,
    guides,
    labs,
    scale_y_continuous,
    theme,
)
from plotnine.stats.stat_sina import stat_sina

from ._aggregate import tidy_expression
from ._compat import SKIP_GGANN_FRAMES, accepts_groupby
from ._expression import ordered_unique, resolve_source, source_var_names
from ._grouping import _downsample_cells, _drop_missing_groups, _group_categories, _order_groups
from ._palette import scale_color_obs, scale_fill_obs
from ._resolve import Ref, parse_token, plain_name, resolve_frame
from .plots import (
    _continuous_scale,
    _feature_facet,
    _is_numeric,
)
from .publication import _family_theme

__all__ = [
    "plot_violin",
    "plot_stacked_violin",
    "plot_box",
    "plot_sina",
    "plot_ridge",
    "plot_expression_bar",
    "plot_expression_line",
]

# A ridge draws two layers per group; ~13 s at 300 groups, so refuse past this.
MAX_RIDGE_GROUPS = 100


def _positive_int(name: str, value, minimum: int = 1) -> int:
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)):
        raise ValueError(f"{name} must be an integer >= {minimum}, got {value!r}.")
    if value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}, got {value!r}.")
    return int(value)


def _positive_finite(name: str, value) -> float:
    ok = isinstance(value, (int, float, np.integer, np.floating))
    if not ok or isinstance(value, (bool, np.bool_)) or not math.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be a finite number > 0, got {value!r}.")
    return float(value)


def _jitter_rng(random_state):
    """A fresh seeded generator per plot (seed 0 included); ``None`` stays nondeterministic.

    plotnine treats an int seed of 0 as "unseeded" and keeps a passed-in generator's
    state only on a deep copy, so hand each layer its own ``RandomState``.
    """
    return None if random_state is None else np.random.RandomState(random_state)


class _stat_sina_seeded(stat_sina):
    """``stat_sina`` that restarts from its seed on every build.

    plotnine replaces the ``random_state`` parameter by a live generator the first
    time it builds, so a second ``draw()``/export of the same plot would continue
    the stream and move the points (and a ``None`` seed leaves the ``numpy.random``
    module there, which cannot be deep-copied by a later ``+``).
    """

    def __init__(self, mapping=None, data=None, **kwargs):
        super().__init__(mapping, data, **kwargs)
        self._seed = self.params["random_state"]

    def setup_params(self, data):
        seed = self._seed
        self.params["random_state"] = np.random.RandomState(seed)  # None: fresh entropy
        super().setup_params(data)


def _select_order(adata, group_by: str, order, categories_order):
    """Subset to the groups in ``order`` (scanpy's ``order=``) and use it as the axis order."""
    if order is None:
        return adata, categories_order
    if categories_order is not None:
        raise ValueError("Pass order= (a subset) or categories_order= (the full order), not both.")
    order = ordered_unique(order)
    present = set(adata.obs[group_by].dropna().astype(object))
    unknown = [g for g in order if g not in present]
    if not order or unknown:
        raise ValueError(
            f"order must be non-empty and name groups of {group_by!r}; unknown: {unknown}."
        )
    return adata[adata.obs[group_by].isin(order).to_numpy()], order


def _log1p_values(tidy: pd.DataFrame) -> pd.DataFrame:
    """``log1p`` the plotted values (zeros stay zero, unlike a log y scale)."""
    if np.any(tidy["value"].to_numpy() <= -1):
        raise ValueError(
            "log=True applies log1p, which needs values > -1; "
            "use unscaled expression (a log-normalised layer or adata.raw)."
        )
    tidy["value"] = np.log1p(tidy["value"])
    return tidy


def _tidy_features(adata, features, group_by, *, layer=None, use_raw=None, extra_obs=()):
    """Long per-cell table for genes and/or numeric obs columns; returns ``(tidy, names)``.

    A bare name is a gene when the expression source has it (so forced-gene
    behaviour matches the matrix helpers), else a numeric ``adata.obs`` column.
    ``obs("x")`` / ``"obs:x"`` forces obs and ``gene("x")`` forces a gene.
    """
    features = ordered_unique(features)
    if not features:
        raise ValueError("No genes requested; pass at least one gene name.")
    kind, _ = resolve_source(adata, layer, use_raw)
    var = source_var_names(adata, kind)
    names, from_obs = [], []
    for item in features:
        tok = parse_token(item)
        if not isinstance(tok, (Ref, str)):
            raise ValueError(f"{item!r} is not a gene or obs column selector.")
        if isinstance(tok, Ref):
            if tok.source == "gene" and (tok.layer is not None or tok.use_raw is not None):
                raise ValueError(
                    f"{tok!r}: per-gene layer/use_raw is not supported here; pass layer=/use_raw=."
                )
            name, is_obs = tok.name, tok.source == "obs"
        else:
            name, is_obs = tok, tok not in var and tok in adata.obs.columns
        names.append(name)
        from_obs.append(is_obs)
    if len(set(names)) < len(names):
        raise ValueError(f"Features resolve to duplicate names: {names}.")
    obs_names = [n for n, o in zip(names, from_obs) if o]
    missing = [n for n in obs_names if n not in adata.obs.columns]
    if missing:
        raise KeyError(f"Observation column(s) not found: {', '.join(map(repr, missing))}.")
    non_numeric = [n for n in obs_names if not _is_numeric(adata.obs[n])]
    if non_numeric:
        raise TypeError(
            f"Observation column(s) {non_numeric} are not numeric; cannot plot a distribution."
        )

    gene_names = [n for n, o in zip(names, from_obs) if not o]
    parts = []
    if gene_names:
        parts.append(
            tidy_expression(
                adata, gene_names, group_by, layer=layer, use_raw=use_raw, extra_obs=extra_obs
            )
        )
    if not obs_names:
        return parts[0], names
    meta = list(dict.fromkeys([group_by, *extra_obs]))
    for n in obs_names:
        frame = adata.obs[meta].reset_index(drop=True)
        frame.insert(0, "obs_name", pd.Series(adata.obs_names))
        frame.insert(1, "feature", n)
        frame.insert(2, "value", adata.obs[n].to_numpy(dtype=float))
        parts.append(frame)
    tidy = pd.concat(parts, ignore_index=True)
    tidy["feature"] = pd.Categorical(tidy["feature"].astype(str), categories=names, ordered=True)
    return tidy, names


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


def _add_stats(plot):
    """Group-comparison label 10 % of the data range above the data, with headroom for it.

    The default label sits at 95 % of the range, centred on the middle group, where
    it overprints a tall violin tip or box whisker.
    """
    return (
        plot
        + pe.stat_compare_means(label_y_npc=1.1)
        + scale_y_continuous(expand=(0.05, 0, 0.18, 0))
    )


@accepts_groupby
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
    jitter_size: float = 0.6,
    jitter_alpha: float = 0.5,
    stats: bool = False,
    categories_order: Sequence[str] | None = None,
    order: Sequence[str] | None = None,
    downsample: int | None = None,
    random_state: int | None = 0,
    log: bool = False,
    rotation: float = 45,
):
    """Per-group expression box plots, one facet per gene, with jittered cells overlaid.

    Set ``jitter=False`` for a plain box plot or ``split_by`` for a gene x split
    facet grid (``ncol`` only applies without ``split_by``). ``genes`` may also
    name numeric ``adata.obs`` columns (see :func:`ggann.plot_violin`). ``stats=True`` adds one two-sided, unpaired global comparison per
    facet: Mann-Whitney U for two groups and Kruskal-Wallis for more than two.
    The label is the three-significant-digit p-value; no multiplicity correction
    is applied across facets. Pass ``downsample=N`` to cap cells per group before
    the (jitter) draw for large data — see :func:`ggann.plot_violin`.

    Note: ``downsample`` subsets the cells the geoms see, so with ``stats=True``
    the p-value is computed on the subsample, and the boxplot's outliers reflect
    only the kept cells. ``random_state`` seeds the downsampling and the point
    jitter, so the same seed (0 included) gives identical point positions on every
    draw and export; ``random_state=None`` is explicitly nondeterministic.

    Parameters
    ----------
    adata
        Annotated data matrix.
    genes : sequence of str
        Genes (or numeric obs columns) shown in facets.
    group_by : str
        Observation grouping column.
    split_by : str, optional
        Additional facet column.
    layer, use_raw : optional
        Mutually exclusive expression source.
    ncol : int
        Facet columns; ignored when ``split_by`` is given.
    jitter, stats : bool
        Add cell points or a group-comparison layer.
    jitter_size, jitter_alpha : float
        Point size and opacity.
    categories_order : sequence of str, optional
        Complete observed group order.
    order : sequence of str, optional
        Show only these groups, in this order (scanpy ``order=``).
    downsample : int, optional
        Maximum cells per group.
    random_state : int, optional
        Downsampling and jitter seed; ``None`` is nondeterministic.
    log : bool
        Plot ``log1p(expression)``; zeros stay at zero.
    rotation : float
        Group-label angle in degrees.

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
    adata = _drop_missing_groups(adata, group_by, split_by)
    adata, categories_order = _select_order(adata, group_by, order, categories_order)
    adata = _downsample_cells(adata, group_by, downsample, random_state=random_state)
    extra = [split_by] if split_by else []
    tidy, genes = _tidy_features(
        adata, genes, group_by, layer=layer, use_raw=use_raw, extra_obs=extra
    )
    if log:
        tidy = _log1p_values(tidy)
    if categories_order is None:
        categories_order = _group_categories(adata, group_by)
    tidy = _order_groups(tidy, group_by, categories_order)

    plot = ggplot(tidy, aes(group_by, "value", fill=group_by)) + geom_boxplot(
        width=0.7, outlier_alpha=0.0 if jitter else 1.0
    )
    if jitter:
        plot = plot + geom_jitter(
            width=0.2,
            height=0.0,
            size=jitter_size,
            alpha=jitter_alpha,
            stroke=0,
            random_state=_jitter_rng(random_state),
        )
    plot = (
        plot
        + _feature_facet(split_by, ncol=ncol)
        + scale_fill_obs(adata, group_by)
        + labs(x="", y="log1p(expression)" if log else "expression", fill=group_by)
        + _family_theme("distribution")
        + pe.rotate_x_text(rotation)
    )
    return _add_stats(plot) if stats else plot


@accepts_groupby
def plot_sina(
    adata,
    genes: Sequence[str],
    group_by: str,
    *,
    split_by: str | None = None,
    layer: str | None = None,
    use_raw: bool | None = None,
    ncol: int = 1,
    size: float = 0.8,
    alpha: float = 0.6,
    violin: bool = True,
    bins: int = 50,
    categories_order: Sequence[str] | None = None,
    order: Sequence[str] | None = None,
    downsample: int | None = None,
    random_state: int | None = 0,
    log: bool = False,
    rotation: float = 45,
):
    """Sina / beeswarm of per-group expression, one facet per gene.

    A sina plot spreads each group's cells horizontally in proportion to the local
    density (via plotnine-extra's ``geom_sina``), so it shows every cell like a
    jitter but keeps the shape of a violin. ``violin=True`` draws a faint violin
    behind the points for context. ``downsample=N`` caps cells per group first --
    useful for large data, since a sina draws one mark per cell. ``random_state``
    seeds the downsampling and the horizontal spread, so a fixed seed (0
    included) repeats the same positions on every draw; ``None`` is nondeterministic.

    Known limitation: plotnine's ``stat_sina`` also jitters the *y* values of a
    panel when every value in it is an integer (for example raw counts), so
    plotted heights can move by up to about one unit. ggann does not alter the
    data to hide this; use :func:`ggann.plot_violin` or :func:`ggann.plot_box`
    for integer-valued data, or a normalised layer. See
    ``plans/evaluation/upstream/U5.md``.

    Parameters
    ----------
    adata
        Annotated data matrix.
    genes : sequence of str
        Genes (or numeric obs columns) shown in facets.
    group_by : str
        Observation grouping column.
    split_by : str, optional
        Additional facet column.
    layer, use_raw : optional
        Mutually exclusive expression source.
    ncol : int
        Facet columns; ignored when ``split_by`` is given.
    size, alpha : float
        Point size and opacity.
    violin : bool
        Draw a contextual violin.
    bins : int
        Density bins used by the sina geometry; a positive integer.
    categories_order : sequence of str, optional
        Complete observed group order.
    order : sequence of str, optional
        Show only these groups, in this order (scanpy ``order=``).
    downsample : int, optional
        Maximum cells per group.
    random_state : int, optional
        Downsampling and sina-spread seed; ``None`` is nondeterministic.
    log : bool
        Plot ``log1p(expression)``; zeros stay at zero.
    rotation : float
        Group-label angle in degrees.

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
    bins = _positive_int("bins", bins)
    adata = _drop_missing_groups(adata, group_by, split_by)
    adata, categories_order = _select_order(adata, group_by, order, categories_order)
    adata = _downsample_cells(adata, group_by, downsample, random_state=random_state)
    extra = [split_by] if split_by else []
    tidy, genes = _tidy_features(
        adata, genes, group_by, layer=layer, use_raw=use_raw, extra_obs=extra
    )
    if log:
        tidy = _log1p_values(tidy)
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
    ranges = pd.Series(
        tidy.groupby("feature", observed=True)["value"].agg(lambda s: s.max() - s.min())
    )
    ranges = ranges.loc[ranges.to_numpy() > 0]
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
        + pe.geom_sina(
            size=size,
            alpha=alpha,
            stroke=0,
            binwidth=binwidth,
            random_state=random_state,
            stat=_stat_sina_seeded,
        )
        + _feature_facet(split_by, ncol=ncol)
        + scale_color_obs(adata, group_by)
        # the points are tiny; draw larger, opaque legend keys
        + guides(color=guide_legend(override_aes={"size": 3, "alpha": 1}))
        + labs(x="", y="log1p(expression)" if log else "expression", color=group_by)
        + _family_theme("distribution")
        + pe.rotate_x_text(rotation)
    )
    if violin:  # a fill scale without a fill aesthetic warns
        plot = plot + scale_fill_obs(adata, group_by)
    return plot


@accepts_groupby
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
    adata = _drop_missing_groups(adata, group_by, split_by)
    extra = [split_by] if split_by else []
    tidy, genes = _tidy_features(
        adata, genes, group_by, layer=layer, use_raw=use_raw, extra_obs=extra
    )
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
        + _family_theme("distribution")
        + pe.rotate_x_text(45)
    )
    if error != "none":
        plot = plot + geom_errorbar(aes(ymin="ymin", ymax="ymax"), width=0.3)
    return plot


@accepts_groupby
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
    bins: int | None = None,
):
    """Aggregated expression trend across an ordered variable ``x``, one facet per gene.

    ``x`` is an obs column (an ordered condition, timepoint or a numeric covariate
    such as pseudotime). With ``group_by`` given, one line is drawn per group and
    coloured by it (e.g. an expression trajectory per cell type across timepoints).
    ``agg`` sets the summarised value (``"mean"`` default, or ``"median"`` etc.);
    ``error`` adds a summary error bar per point (``"se"`` / ``"sd"`` / ``"none"``).
    A continuous ``x`` (pseudotime) summarises every distinct value, which is one
    cell per point; pass ``bins=N`` to cut it into ``N`` equal-width bins plotted
    at the bin centres.

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
    bins : int, optional
        Number of equal-width bins for a numeric ``x``.

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
    if bins is not None:
        bins = _positive_int("bins", bins)
    adata = _drop_missing_groups(adata, group_by)
    genes = ordered_unique(genes)
    if not genes:
        raise ValueError("No genes requested; pass at least one gene name.")
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
    frame = resolve_frame(
        adata, cols, layer=layer, use_raw=use_raw, strict=True
    )  # already densified
    long = frame.melt(
        id_vars=id_vars, value_vars=gene_names, var_name="feature", value_name="value"
    )
    long["feature"] = pd.Categorical(long["feature"], categories=gene_names, ordered=True)

    if bins is not None:
        if not _is_numeric(long[xname]):
            raise ValueError(f"bins needs a numeric x; {xname!r} is not numeric.")
        cut, edges = pd.cut(long[xname], bins, retbins=True)
        codes = cut.cat.codes.to_numpy()
        edges = np.asarray(edges, dtype=float)
        centres = (edges[:-1] + edges[1:]) / 2
        long[xname] = np.where(codes >= 0, centres[codes], np.nan)

    # Order the x axis: numeric stays numeric; categorical keeps its obs order.
    if not _is_numeric(long[xname]):
        x_cats = _group_categories(adata, x)
        if x_cats is None:
            x_cats = list(pd.unique(pd.Series(long[xname]).to_numpy()))
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
        + _family_theme("distribution")
    )
    if color is not None:
        plot = plot + scale_color_obs(adata, gname)
    if error != "none":
        plot = plot + geom_errorbar(aes(ymin="ymin", ymax="ymax"), width=0.15)
    return plot


@accepts_groupby
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
    order: Sequence[str] | None = None,
    log: bool = False,
    rotation: float = 45,
    jitter_size: float = 0.8,
    jitter_alpha: float = 0.5,
):
    """Per-group expression distributions, one facet per gene (stacked-violin style).

    ``add_box=True`` (default) nests a slim white boxplot inside each violin so the
    median and quartiles read off cleanly, the way scplotter's ``FeatureStatPlot``
    does. ``add_points=True`` overlays the individual cells as jitter (scplotter's
    ``add_point``; scanpy's ``stripplot``). ``split_by`` adds a facet column (gene
    rows x split columns; ``ncol`` then has no effect). Entries of ``genes`` may
    also be numeric ``adata.obs`` columns such as ``"n_genes"`` (a bare name is a
    gene when the expression source has it; use ``ggann.obs("n_genes")`` to force
    the obs column, as ``sc.pl.violin`` does with ``keys``).
    Set ``stats=True`` to add one two-sided, unpaired global comparison per
    facet: Mann-Whitney U for two groups and Kruskal-Wallis for more than two.
    The label is the three-significant-digit p-value; no multiplicity correction
    is applied across facets. ``downsample`` caps cells per group before the
    (slow) violin KDE -- a big speed-up on large data for a visually identical
    plot.

    Note: ``downsample`` subsets the cells the geoms see, so any ``stats=True``
    p-value is then computed on the *subsample*, not the full data. Leave
    ``downsample`` unset when you need the reported test to reflect every cell.
    ``random_state`` controls which cells are retained and seeds the point jitter:
    a fixed seed (0 included) gives identical point positions on every draw and
    export, ``None`` is nondeterministic.

    Parameters
    ----------
    adata
        Annotated data matrix.
    genes : sequence of str
        Genes (or numeric obs columns) shown in facets.
    group_by : str
        Observation grouping column.
    split_by : str, optional
        Additional facet column.
    layer, use_raw : optional
        Mutually exclusive expression source.
    ncol : int
        Facet columns; ignored when ``split_by`` is given.
    scale : str
        plotnine violin-width scaling mode.
    add_box, add_points, stats : bool
        Add inner boxes, cell points, or group tests.
    downsample : int, optional
        Maximum cells per group.
    random_state : int, optional
        Downsampling and jitter seed; ``None`` is nondeterministic.
    categories_order : iterable of str, optional
        Complete observed group order.
    order : sequence of str, optional
        Show only these groups, in this order (scanpy ``order=``).
    log : bool
        Plot ``log1p(expression)``; zeros stay at zero (never a log axis).
    rotation : float
        Group-label angle in degrees.
    jitter_size, jitter_alpha : float
        Size and opacity of the ``add_points`` cells.

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
    adata = _drop_missing_groups(adata, group_by, split_by)
    adata, categories_order = _select_order(adata, group_by, order, categories_order)
    adata = _downsample_cells(adata, group_by, downsample, random_state=random_state)
    extra = [split_by] if split_by else []
    tidy, genes = _tidy_features(
        adata, genes, group_by, layer=layer, use_raw=use_raw, extra_obs=extra
    )
    if log:
        tidy = _log1p_values(tidy)
    if categories_order is None:
        categories_order = _group_categories(adata, group_by)
    tidy = _order_groups(tidy, group_by, categories_order)
    plot = ggplot(tidy, aes(group_by, "value", fill=group_by)) + geom_violin(scale=scale)
    # box first, then points on top -- otherwise the white box occludes the jitter
    if add_box:
        plot = plot + geom_boxplot(width=0.12, fill="white", outlier_alpha=0.0, show_legend=False)
    if add_points:
        plot = plot + geom_jitter(
            width=0.2,
            height=0.0,
            size=jitter_size,
            alpha=jitter_alpha,
            stroke=0,
            random_state=_jitter_rng(random_state),
        )
    plot = (
        plot
        + _feature_facet(split_by, ncol=ncol)
        + scale_fill_obs(adata, group_by)
        + labs(x="", y="log1p(expression)" if log else "expression", fill=group_by)
        + _family_theme("distribution")
        + pe.rotate_x_text(rotation)
    )
    return _add_stats(plot) if stats else plot


@accepts_groupby
def plot_stacked_violin(
    adata,
    genes: Sequence[str],
    group_by: str,
    *,
    layer: str | None = None,
    use_raw: bool | None = None,
    scale: str = "width",
    categories_order=None,
    order: Sequence[str] | None = None,
    downsample: int | None = None,
    random_state: int | None = 0,
    log: bool = False,
    rotation: float = 90,
    stripplot: bool = False,
    yticklabels: bool = False,
    cmap: str = "Blues",
):
    """Compact genes-as-rows violin grid across groups (``sc.pl.stacked_violin``).

    Genes are rows and groups run along x, the layout of
    ``sc.pl.stacked_violin(..., swap_axes=True)`` (Scanpy's default puts groups
    as rows). As in Scanpy, each violin is filled by that group's **median**
    expression of the gene (``cmap``, shared colour bar). Each row has its own y range and no y tick labels (``yticklabels``
    restores them), so many genes stay readable. ``downsample`` caps cells per
    group before the KDE; ``random_state`` controls which cells are retained and
    seeds the ``stripplot`` jitter.

    Parameters
    ----------
    adata
        Annotated data matrix.
    genes : sequence of str
        Genes (or numeric obs columns) shown as rows.
    group_by : str
        Observation grouping column.
    layer, use_raw : optional
        Mutually exclusive expression source.
    scale : str
        plotnine violin-width scaling mode.
    categories_order : iterable of str, optional
        Complete observed group order.
    order : sequence of str, optional
        Show only these groups, in this order (scanpy ``order=``).
    downsample : int, optional
        Maximum cells per group.
    random_state : int, optional
        Downsampling and jitter seed; ``None`` is nondeterministic.
    log : bool
        Plot ``log1p(expression)``; zeros stay at zero.
    rotation : float, default=90
        Group-label angle in degrees. Vertical labels, as in Scanpy, never
        hang past the left edge; at smaller angles a long first label can.
    stripplot : bool
        Overlay the individual cells as jitter.
    yticklabels : bool
        Show per-row y tick labels.
    cmap : str
        Colormap for the median fill.

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
    adata = _drop_missing_groups(adata, group_by)
    adata, categories_order = _select_order(adata, group_by, order, categories_order)
    adata = _downsample_cells(adata, group_by, downsample, random_state=random_state)
    tidy, genes = _tidy_features(adata, genes, group_by, layer=layer, use_raw=use_raw)
    if log:
        tidy = _log1p_values(tidy)
    if categories_order is None:
        categories_order = _group_categories(adata, group_by)
    tidy = _order_groups(tidy, group_by, categories_order)
    tidy["median"] = tidy.groupby(["feature", group_by], observed=True)["value"].transform("median")
    plot = (
        ggplot(tidy, aes(group_by, "value", fill="median"))
        + geom_violin(scale=scale)
        + facet_grid("feature ~ .", scales="free_y")
    )
    if stripplot:
        plot = plot + geom_jitter(
            width=0.2,
            height=0.0,
            size=0.5,
            alpha=0.4,
            stroke=0,
            color="black",
            random_state=_jitter_rng(random_state),
        )
    plot = (
        plot
        + _continuous_scale("fill", tidy["median"], cmap, signed=False)
        + labs(x="", y="", fill="median\n" + ("log1p(expr)" if log else "expression"))
        + _family_theme("distribution")
        + theme(strip_text_y=element_text(angle=0))
        + pe.rotate_x_text(rotation)
    )
    if not yticklabels:
        plot = plot + theme(axis_text_y=element_blank(), axis_ticks_major_y=element_blank())
    return plot


def _kde_curve(values: np.ndarray, grid: np.ndarray, bw_adjust: float = 1.0) -> np.ndarray:
    """KDE density on ``grid``; zeros for degenerate (constant / tiny) groups."""
    from scipy.stats import gaussian_kde

    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size < 2 or np.ptp(values) == 0:
        return np.zeros_like(grid)
    try:
        kde = gaussian_kde(values)
        kde.set_bandwidth(float(kde.factor) * bw_adjust)
        dens = kde(grid)
    except np.linalg.LinAlgError:  # singular covariance
        return np.zeros_like(grid)
    peak = dens.max()
    return dens / peak if peak > 0 else dens


@accepts_groupby
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
    bw_adjust: float = 1.0,
):
    """Ridgeline plot: one density ridge per group, stacked and overlapping, per gene.

    ``scale`` sets ridge height relative to row spacing; values above one overlap
    neighbouring ridges. Groups with fewer than two cells or zero
    variance draw a flat baseline. plotnine has no ridgeline geom, so this builds one
    from a per-group gaussian KDE offset vertically and drawn as ``geom_ribbon``.
    ``bw_adjust`` multiplies the KDE bandwidth (Scott's rule); values below one
    follow the data more closely. A gene whose values are all equal draws flat
    baselines and a warning names it. At most ``MAX_RIDGE_GROUPS`` (100) groups
    are drawn, since each adds two layers; a continuous column used as
    ``group_by`` raises instead of taking minutes to build.

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
        KDE evaluation-grid size; an integer of at least 2.
    categories_order : sequence of str, optional
        Complete observed group order.
    bw_adjust : float
        Bandwidth multiplier; a finite number above zero.

    Returns
    -------
    plotnine.ggplot
        Composable ridgeline plot.

    Raises
    ------
    KeyError
        If genes, grouping metadata, or a layer are missing.
    ValueError
        If source, ordering, or KDE configuration is invalid, ``scale`` is not a
        finite positive number, or ``group_by`` has more than
        ``MAX_RIDGE_GROUPS`` levels.

    Notes
    -----
    KDE work scales with genes, groups, cells, and ``n_grid``. Input is unchanged.

    Examples
    --------
    >>> p = plot_ridge(adata, ["CD3D"], group_by="cell_type")
    """
    n_grid = _positive_int("n_grid", n_grid, minimum=2)
    scale = _positive_finite("scale", scale)
    bw_adjust = _positive_finite("bw_adjust", bw_adjust)
    adata = _drop_missing_groups(adata, group_by)
    if group_by in adata.obs:
        n_levels = int(adata.obs[group_by].nunique())
        if n_levels > MAX_RIDGE_GROUPS:
            raise ValueError(
                f"group_by={group_by!r} has {n_levels} levels; plot_ridge draws at most "
                f"{MAX_RIDGE_GROUPS} (two layers per group). Bin or subset the field first."
            )
    tidy, genes = _tidy_features(adata, genes, group_by, layer=layer, use_raw=use_raw)
    if categories_order is None:
        categories_order = _group_categories(adata, group_by)
    tidy = _order_groups(tidy, group_by, categories_order)
    order = list(tidy[group_by].cat.categories)
    pos = {g: i for i, g in enumerate(order)}

    flat = [
        str(f) for f, v in tidy.groupby("feature", observed=True)["value"] if int(v.nunique()) < 2
    ]
    if flat:
        warnings.warn(
            f"plot_ridge: constant values for {flat}; their panels draw flat baselines.",
            UserWarning,
            skip_file_prefixes=SKIP_GGANN_FRAMES,
        )
    rows = []
    for feature, fdf in tidy.groupby("feature", observed=True):
        lo, hi = float(fdf["value"].min()), float(fdf["value"].max())
        if hi <= lo:
            hi = lo + 1.0
        grid = np.linspace(lo, hi, n_grid)
        for g, gdf in fdf.groupby(group_by, observed=True):
            dens = _kde_curve(gdf["value"].to_numpy(), grid, bw_adjust) * scale
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
        + _family_theme("distribution")
        + theme(panel_grid=element_blank())
    )
