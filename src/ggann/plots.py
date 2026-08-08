"""High-level plotnine-native plotting helpers mirroring scanpy's core figures.

Each helper follows the same contract: extract a tidy DataFrame with annplyr,
then hand it to plotnine / plotnine-extra. None of them index ``adata`` directly.
All return ordinary :class:`plotnine.ggplot` objects, so they remain fully
composable with ``+ scale_*`` / ``+ theme(...)`` / ``+ facet_*``.
"""

from __future__ import annotations

import warnings
from typing import Iterable, Sequence

import pandas as pd
import plotnine_extra as pe
from plotnine import (
    aes,
    element_blank,
    element_text,
    facet_grid,
    geom_point,
    geom_tile,
    ggplot,
    guide_legend,
    guides,
    labs,
    scale_color_cmap,
    scale_fill_cmap,
    scale_size,
    theme,
)

from ._aggregate import aggregate_expression, aggregate_means, tidy_expression
from ._expression import ordered_unique
from ._grouping import _downsample_cells, _group_categories, _order_groups
from ._palette import scale_color_obs
from ._resolve import embedding_coords, plain_name, resolve_frame
from .theme import theme_ggann

__all__ = [
    "plot_embedding",
    "plot_features",
    "plot_dotplot",
    "plot_matrixplot",
    "plot_embedding_density",
    "plot_heatmap",
]


def _is_numeric(series: pd.Series) -> bool:
    return pd.api.types.is_numeric_dtype(series) and not isinstance(
        series.dtype, pd.CategoricalDtype
    )


def _feature_facet(split_by: str | None, *, ncol: int = 1, scales: str = "free_y"):
    """Facet by feature, adding a ``split_by`` column dimension when given.

    ``split_by=None`` -> one wrapped panel per gene; ``split_by`` set -> a grid of
    gene (rows) x split level (columns), like scplotter's ``FeatureStatPlot(split_by=)``.
    """
    if split_by is None:
        return pe.facet_wrap("~feature", ncol=ncol, scales=scales)
    return pe.facet_grid(f"feature ~ {split_by}", scales=scales)


def _embedding_axes():
    """Hide the tick numbers on embeddings -- UMAP/t-SNE units are arbitrary,

    so the numbers are noise (scplotter's ``CellDimPlot`` / Seurat drop them too).
    """
    return theme(axis_text=element_blank(), axis_ticks=element_blank())


def _centroid_labels(
    df: pd.DataFrame, cname: str, xcol: str, ycol: str, split_by: str | None = None
) -> pd.DataFrame:
    """Median position of each category, for placing a cluster label at its centre.

    When ``split_by`` is given the centroids are computed *within* each facet, so a
    label lands where its category actually sits in that panel rather than at the
    pooled median (which would be broadcast to every facet).
    """
    keys = [cname] if split_by is None else [split_by, cname]
    cents = df.groupby(keys, observed=True)[[xcol, ycol]].median().reset_index()
    return cents.rename(columns={cname: "label"})


def plot_embedding(
    adata,
    basis: str = "umap",
    color: str | None = None,
    *,
    split_by: str | None = None,
    layer: str | None = None,
    use_raw: bool | None = None,
    size: float = 1.5,
    alpha: float = 0.9,
    pointdensity: bool | None = None,
    label: bool = False,
    label_size: float = 9,
    low: str = "#d9d9d9",
    high: str = "#2166ac",
    downsample: int | None = None,
    random_state: int | None = 0,
):
    """Scatter over an embedding (UMAP/t-SNE/PCA), optionally coloured and split.

    * ``color=None`` -> density-coloured scatter (``geom_pointdensity``), which
      reads well for dense embeddings.
    * categorical ``color`` (e.g. an obs cluster column) -> discrete colours
      (reusing scanpy's stored palette when present).
    * numeric ``color`` (a gene or continuous obs column) -> a gradient.
    * ``split_by`` -> facet the scatter over an obs column (Seurat ``split.by``).
    * ``label=True`` -> for a categorical ``color``, print each category at its
      centroid using repelled (non-overlapping) text, like scplotter's
      ``CellDimPlot`` / Seurat ``label=TRUE``.
    * ``downsample=N`` -> randomly keep at most ``N`` cells before drawing, for a
      much lighter scatter on large data (density reads the same; exact points differ).
      ``random_state`` controls which cells are retained.
    """
    adata = _downsample_cells(adata, None, downsample, random_state=random_state)
    coords = embedding_coords(adata, basis)
    if coords.shape[1] < 2:
        raise ValueError(
            f"Embedding '{basis}' has only {coords.shape[1]} dimension(s); "
            "plot_embedding requires at least 2."
        )
    xcol, ycol = coords.columns[:2]

    split_col = None
    if split_by is not None:
        split_col = resolve_frame(adata, [split_by])[[split_by]]

    def _facet(plot):
        return plot + pe.facet_wrap("~" + split_by) if split_by is not None else plot

    if label and color is None:
        warnings.warn(
            "plot_embedding: label=True is ignored when color is None "
            "(centroid labels need a categorical color).",
            stacklevel=2,
        )

    if color is None:
        if pointdensity is None:
            pointdensity = True
        base = coords if split_col is None else coords.join(split_col)
        if pointdensity:
            return _facet(
                ggplot(base, aes(xcol, ycol))
                + pe.geom_pointdensity(size=size, alpha=alpha)
                + labs(color="density")
                + theme_ggann()
                + _embedding_axes()
            )
        return _facet(
            pe.DimPlot(base, x=xcol, y=ycol, size=size, alpha=alpha)
            + theme_ggann()
            + _embedding_axes()
        )

    # `color` may be a bare name, a prefix string ("gene:CD3D@logcounts") or an accessor
    cname = plain_name(adata, color)
    values = resolve_frame(adata, [color], layer=layer, use_raw=use_raw)
    if cname not in values.columns:
        raise KeyError(f"Could not resolve color={color!r} from obs, genes or obsm.")
    df = coords.join(values[[cname]])
    if split_col is not None and split_by not in df.columns:
        df = df.join(split_col)

    if _is_numeric(df[cname]):
        if label:
            warnings.warn(
                f"plot_embedding: label=True is ignored for the numeric color {color!r} "
                "(centroid labels need a categorical color).",
                stacklevel=2,
            )
        # Draw low-expression cells first so high-expression cells are not occluded
        # (mirrors scanpy's sc.pl.embedding ordering).
        df = df.sort_values(cname)
        return _facet(
            pe.FeatureDimPlot(
                df,
                feature=cname,
                x=xcol,
                y=ycol,
                low=low,
                high=high,
                size=size,
                alpha=alpha,
            )
            + theme_ggann()
            + _embedding_axes()
        )
    plot = (
        pe.DimPlot(df, x=xcol, y=ycol, color=cname, size=size, alpha=alpha)
        + scale_color_obs(adata, cname)
        # enlarge the legend swatches so categories stay readable (scplotter does this)
        + guides(color=guide_legend(override_aes={"size": 4}))
        + theme_ggann()
        + _embedding_axes()
    )
    if label:
        cents = _centroid_labels(df, cname, xcol, ycol, split_by=split_by)
        # white-backed repelled labels at centroids, like scplotter's label_bg="white"
        plot = plot + pe.geom_label_repel(
            aes(xcol, ycol, label="label"),
            data=cents,
            size=label_size * 0.85,
            fill="white",
            color="black",
            inherit_aes=False,
        )
    return _facet(plot)


def plot_features(
    adata,
    features: Sequence[str],
    basis: str = "umap",
    *,
    layer: str | None = None,
    use_raw: bool | None = None,
    ncol: int | None = None,
    size: float = 1.2,
    alpha: float = 0.9,
    cmap: str = "magma",
    downsample: int | None = None,
    random_state: int | None = 0,
):
    """Multi-gene embedding grid: one faceted panel per feature.

    Panels share a single expression colour scale, which is best for comparing
    magnitudes *across* genes. Because the scale is shared, a low-range gene shown
    next to a high-range one will look faint -- for independent per-gene colour
    bars (like ``sc.pl.umap(color=[...])``), compose separate ``plot_embedding``
    calls with the re-exported ``Wrap`` / ``plot_layout`` instead. ``downsample=N``
    caps cells before drawing for lighter panels on large data; ``random_state``
    controls the retained cells.
    """
    features = ordered_unique(features)
    adata = _downsample_cells(adata, None, downsample, random_state=random_state)
    coords = embedding_coords(adata, basis)
    if coords.shape[1] < 2:
        raise ValueError(f"Embedding '{basis}' has fewer than 2 dimensions.")
    xcol, ycol = coords.columns[:2]

    values = resolve_frame(adata, list(features), layer=layer, use_raw=use_raw)
    feats = list(
        dict.fromkeys(
            f for f in features if f in values.columns and _is_numeric(values[f])
        )
    )
    if not feats:
        raise ValueError(
            "plot_features needs at least one numeric feature (gene or metric)."
        )

    df = coords.join(values[feats])
    long = df.melt(
        id_vars=[xcol, ycol],
        value_vars=feats,
        var_name="feature",
        value_name="expression",
    ).sort_values("expression")
    long["feature"] = pd.Categorical(long["feature"], categories=feats, ordered=True)
    return (
        ggplot(long, aes(xcol, ycol, color="expression"))
        + geom_point(size=size, alpha=alpha)
        + pe.facet_wrap("~feature", ncol=ncol)
        + scale_color_cmap(cmap_name=cmap)
        + theme_ggann()
        + _embedding_axes()
    )


def _cell_rank(tidy: pd.DataFrame, group_by: str) -> pd.DataFrame:
    """Add a ``cell_rank`` column that orders cells by their ``group_by`` category.

    Cells are sorted by group then numbered 0..N-1, so a per-cell x layout reads
    left-to-right by group. Shared by :func:`plot_heatmap` and
    ``markers.plot_tracksplot``.
    """
    positions = tidy.groupby("feature", observed=True, sort=False).cumcount()
    cell = tidy[["obs_name", group_by]].drop_duplicates()
    n_cells = int(positions.max()) + 1 if len(positions) else 0
    if cell["obs_name"].is_unique and len(cell) == n_cells:
        cell = cell.sort_values(group_by)
        ranks = pd.Series(range(len(cell)), index=cell["obs_name"])
        result = tidy.copy()
        result["cell_rank"] = result["obs_name"].map(ranks).to_numpy()
        return result

    first_feature = ~positions.duplicated()
    cell = pd.DataFrame(
        {
            "position": positions[first_feature].to_numpy(),
            group_by: tidy.loc[first_feature, group_by].to_numpy(),
        }
    ).sort_values(group_by, kind="stable")
    ranks = pd.Series(range(len(cell)), index=cell["position"])
    result = tidy.copy()
    result["cell_rank"] = positions.map(ranks).to_numpy()
    return result


def plot_dotplot(
    adata,
    genes: Sequence[str],
    group_by: str,
    *,
    split_by: str | None = None,
    layer: str | None = None,
    use_raw: bool | None = None,
    standard_scale: str | None = None,
    expression_cutoff: float = 0.0,
    cmap: str = "Reds",
    size_range: tuple[float, float] = (0.5, 8.0),
    categories_order: Iterable[str] | None = None,
):
    """Marker dotplot: dot *size* = fraction expressing, *colour* = mean expression.

    Defaults to ``adata.raw`` so values match ``sc.pl.dotplot``. ``split_by`` adds a
    facet column so the dotplot is repeated per split level (scplotter ``split_by``).
    """
    genes = ordered_unique(genes)
    extra = [split_by] if split_by else []
    agg = aggregate_expression(
        adata,
        genes,
        group_by,
        layer=layer,
        use_raw=use_raw,
        expression_cutoff=expression_cutoff,
        standard_scale=standard_scale,
        extra_by=extra,
    )
    if categories_order is None:
        categories_order = _group_categories(adata, group_by)
    agg = _order_groups(agg, group_by, categories_order)
    color_lab = "scaled\nexpression" if standard_scale else "mean\nexpression"
    plot = (
        ggplot(agg, aes("feature", group_by))
        + geom_point(aes(size="fraction", color="mean_expression"))
        + scale_color_cmap(cmap_name=cmap)
        + scale_size(range=size_range, labels=lambda xs: [f"{x:.0%}" for x in xs])
        + labs(x="", y="", color=color_lab, size="fraction\nexpressing")
        + theme_ggann()
        + pe.rotate_x_text(45)
    )
    return plot + pe.facet_wrap(f"~{split_by}") if split_by else plot


def plot_matrixplot(
    adata,
    genes: Sequence[str],
    group_by: str,
    *,
    split_by: str | None = None,
    layer: str | None = None,
    use_raw: bool | None = None,
    standard_scale: str | None = None,
    cmap: str = "viridis",
    categories_order: Iterable[str] | None = None,
):
    """Aggregated mean-expression heatmap (genes x groups) as a plotnine tile plot.

    Like ``sc.pl.matrixplot``, this defaults to ``standard_scale=None`` (raw group
    means). Pass ``standard_scale="var"`` to rescale each gene to ``[0, 1]``, which
    keeps cross-gene patterns legible when magnitudes differ widely. ``split_by``
    adds a facet column so the heatmap is repeated per split level.
    """
    genes = ordered_unique(genes)
    extra = [split_by] if split_by else []
    agg = aggregate_means(
        adata,
        genes,
        group_by,
        layer=layer,
        use_raw=use_raw,
        standard_scale=standard_scale,
        extra_by=extra,
    )
    if categories_order is None:
        categories_order = _group_categories(adata, group_by)
    agg = _order_groups(agg, group_by, categories_order)
    color_lab = "scaled\nexpression" if standard_scale else "mean\nexpression"
    plot = (
        ggplot(agg, aes("feature", group_by, fill="mean_expression"))
        + geom_tile()
        + scale_fill_cmap(cmap_name=cmap)
        + labs(x="", y="", fill=color_lab)
        + theme_ggann()
        + pe.rotate_x_text(45)
    )
    return plot + pe.facet_wrap(f"~{split_by}") if split_by else plot


def plot_embedding_density(
    adata,
    basis: str = "umap",
    group_by: str | None = None,
    *,
    size: float = 2.0,
    ncol: int | None = None,
    cmap: str = "viridis",
    downsample: int | None = None,
    random_state: int | None = 0,
):
    """Per-group cell density over an embedding.

    For each ``group_by`` category a 2D Gaussian KDE is fit on that group's
    embedding coordinates and evaluated at its cells, so every panel shows *where
    that group's cells concentrate* on the shared embedding, min-max scaled to
    ``[0, 1]``. With ``group_by=None`` a single density over all cells is drawn.

    This computes the density directly (via ``scipy.stats.gaussian_kde``) rather
    than reading a pre-computed ``sc.tl.embedding_density`` result, so it is a
    native alternative rather than a byte-for-byte reproduction of scanpy's output.
    ``downsample`` caps cells per group; ``random_state`` controls the sample.
    """
    import numpy as np
    from scipy.stats import gaussian_kde

    adata = _downsample_cells(adata, group_by, downsample, random_state=random_state)
    coords = embedding_coords(adata, basis)
    if coords.shape[1] < 2:
        raise ValueError(f"Embedding '{basis}' has fewer than 2 dimensions.")
    xcol, ycol = coords.columns[:2]

    def _density(sub: pd.DataFrame) -> pd.Series:
        xy = sub[[xcol, ycol]].to_numpy().T  # shape (2, n_cells)
        # KDE needs >2 points and some spread; degenerate groups get a flat density
        if xy.shape[1] < 3 or float(xy.std(axis=1).min()) == 0.0:
            return pd.Series(0.0, index=sub.index)
        try:
            d = gaussian_kde(xy)(xy)
        except np.linalg.LinAlgError:  # singular covariance -- fall back to flat
            return pd.Series(0.0, index=sub.index)
        lo, hi = float(d.min()), float(d.max())
        d = (d - lo) / (hi - lo) if hi > lo else d * 0.0
        return pd.Series(d, index=sub.index)

    if group_by is None:
        df = coords.copy()
        df["density"] = _density(df).to_numpy()
        plot = ggplot(df, aes(xcol, ycol, color="density")) + geom_point(size=size)
    else:
        gcol = resolve_frame(adata, [group_by])[[group_by]]
        df = coords.join(gcol)
        # groupby(...).apply concatenates in group-sorted order, so reindex back
        # to df's cell order (a bare .to_numpy() would misalign interleaved groups)
        per_group = df.groupby(group_by, observed=True, group_keys=False).apply(
            _density
        )
        df["density"] = per_group.reindex(df.index).to_numpy()
        cats = _group_categories(adata, group_by)
        df = _order_groups(df, group_by, cats)
        plot = (
            ggplot(df, aes(xcol, ycol, color="density"))
            + geom_point(size=size)
            + pe.facet_wrap("~" + group_by, ncol=ncol)
        )
    return plot + scale_color_cmap(cmap_name=cmap) + theme_ggann() + _embedding_axes()


def plot_heatmap(
    adata,
    genes: Sequence[str],
    group_by: str,
    *,
    layer: str | None = None,
    use_raw: bool | None = None,
    cmap: str = "viridis",
    standard_scale: str | None = None,
    categories_order: Sequence[str] | None = None,
    downsample: int | None = None,
    random_state: int | None = 0,
):
    """Per-**cell** expression heatmap, cells grouped along x (``sc.pl.heatmap``).

    One column per cell (blocked and labelled by ``group_by``), one row per gene,
    tile-coloured by expression -- the per-cell counterpart to the aggregated
    :func:`plot_matrixplot`. ``standard_scale='var'`` z-min-max-scales each gene to
    ``[0, 1]`` across cells so low- and high-range genes stay comparable.
    ``downsample=N`` caps cells per group first; ``random_state`` controls which
    cells are retained.
    """
    adata = _downsample_cells(adata, group_by, downsample, random_state=random_state)
    genes = ordered_unique(genes)
    tidy = tidy_expression(adata, genes, group_by, layer=layer, use_raw=use_raw)
    if categories_order is None:
        categories_order = _group_categories(adata, group_by)
    tidy = _order_groups(tidy, group_by, categories_order)

    if standard_scale not in {None, "var", "group", "zscore"}:
        raise ValueError(
            "standard_scale must be None, 'var', 'group' or 'zscore', "
            f"got {standard_scale!r}"
        )
    if standard_scale is not None:
        key = "obs_name" if standard_scale == "group" else "feature"
        g = tidy.groupby(key, observed=True)["value"]
        if standard_scale == "zscore":
            centre = g.transform("mean")
            spread = g.transform(lambda values: values.std(ddof=0)).replace(0, 1)
            tidy["value"] = (tidy["value"] - centre) / spread
        else:
            lo = g.transform("min")
            rng = (g.transform("max") - lo).replace(0, 1)
            tidy["value"] = (tidy["value"] - lo) / rng
        fill_lab = "scaled expr."
    else:
        fill_lab = "expression"

    # order cells by group, then give each a rank so tiles sit side by side
    tidy = _cell_rank(tidy, group_by)
    tidy["feature"] = pd.Categorical(
        tidy["feature"], categories=list(reversed(genes)), ordered=True
    )

    return (
        ggplot(tidy, aes("cell_rank", "feature", fill="value"))
        + geom_tile()
        + facet_grid(". ~ " + group_by, scales="free_x", space="free_x")
        + scale_fill_cmap(cmap_name=cmap)
        + labs(x="", y="", fill=fill_lab)
        + theme_ggann()
        + theme(
            axis_text_x=element_blank(),
            axis_ticks_major_x=element_blank(),
            strip_text_x=element_text(angle=90),
        )
    )
