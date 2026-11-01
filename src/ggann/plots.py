"""High-level plotnine-native plotting helpers mirroring scanpy's core figures.

Each helper follows the same contract: extract a tidy DataFrame with annplyr,
then hand it to plotnine / plotnine-extra. None of them index ``adata`` directly.
All return ordinary :class:`plotnine.ggplot` objects, so they remain fully
composable with ``+ scale_*`` / ``+ theme(...)`` / ``+ facet_*``.
"""

from __future__ import annotations

import warnings
from typing import Iterable, Sequence, cast

import numpy as np
import pandas as pd
import plotnine_extra as pe
from plotnine import (
    aes,
    coord_equal,
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
    scale_color_gradient,
    scale_color_gradient2,
    scale_color_manual,
    scale_fill_cmap,
    scale_fill_gradient2,
    scale_size,
    scale_x_discrete,
    scale_y_discrete,
    scale_y_reverse,
    theme,
)
from scipy import sparse

from ._aggregate import (
    aggregate_expression,
    aggregate_expression_native,
    aggregate_means,
    aggregate_means_native,
    tidy_expression,
)
from ._annotation import annotation_threshold, geom_contrast_text
from ._expression import expression_matrix, ordered_unique, resolve_source, source_var_names
from ._grouping import (
    _check_group_cardinality,
    _downsample_cells,
    _drop_missing_groups,
    _group_categories,
    _order_groups,
)
from ._matplotlib_backend import categorical_palette, promote_matplotlib_plot
from ._palette import obs_colors, scale_color_obs
from ._resolve import (
    ObsmRef,
    Ref,
    embedding_key,
    obsm,
    parse_token,
    plain_name,
    resolve_frame,
)
from .publication import _active_style, _family_theme

__all__ = [
    "plot_embedding",
    "plot_features",
    "plot_dotplot",
    "plot_matrixplot",
    "plot_embedding_density",
    "plot_heatmap",
]


def _is_numeric(series: pd.Series) -> bool:
    """True for continuous values; bool and categorical are discrete."""
    return (
        pd.api.types.is_numeric_dtype(series)
        and not pd.api.types.is_bool_dtype(series)
        and not isinstance(series.dtype, pd.CategoricalDtype)
    )


def _validate_backend(backend: str) -> None:
    if backend not in {"plotnine", "matplotlib"}:
        raise ValueError("backend must be 'plotnine' or 'matplotlib'.")


def _continuous_scale(
    aesthetic: str,
    values: pd.Series,
    cmap: str,
    *,
    signed: bool,
):
    """Publication-aware continuous scale while freezing legacy defaults."""
    style = _active_style()
    scale_cmap = scale_color_cmap if aesthetic == "color" else scale_fill_cmap
    scale_diverging = scale_color_gradient2 if aesthetic == "color" else scale_fill_gradient2
    if style is not None and signed:
        finite = values.to_numpy(dtype=float)
        finite = finite[np.isfinite(finite)]
        limit = float(np.abs(finite).max()) if len(finite) else 1.0
        limit = limit or 1.0
        return scale_diverging(
            low=style.diverging[0],
            mid=style.diverging[1],
            high=style.diverging[2],
            midpoint=0,
            limits=(-limit, limit),
            na_value=style.missing_color,
        )
    if style is not None:
        if cmap == "Reds":
            cmap = style.sequential_cmap
        return scale_cmap(cmap_name=cmap, na_value=style.missing_color)
    return scale_cmap(cmap_name=cmap)


_MAX_DENSITY_GROUPS = 100  # one facet per group; more panels are unreadable and slow


def _kde_density(xy: np.ndarray) -> np.ndarray:
    """Gaussian KDE of ``xy`` (shape ``(2, n)``) at its own points, min-max to [0, 1].

    Degenerate groups (fewer than 3 cells, no spread, singular covariance) get a
    flat density of 0.
    """
    from scipy.stats import gaussian_kde

    if xy.shape[1] < 3 or float(xy.std(axis=1).min()) == 0.0:
        return np.zeros(xy.shape[1])
    try:
        d = gaussian_kde(xy)(xy)
    except np.linalg.LinAlgError:
        return np.zeros(xy.shape[1])
    lo, hi = float(d.min()), float(d.max())
    return (d - lo) / (hi - lo) if hi > lo else np.zeros_like(d)


def _native_vector(matrix) -> np.ndarray:
    if sparse.issparse(matrix):
        return matrix.toarray().reshape(-1)
    return np.asarray(matrix).reshape(-1)


def _native_embedding_frame(
    adata,
    key: str,
    color,
    *,
    layer: str | None,
    use_raw: bool | None,
) -> tuple[pd.DataFrame, str, str, str | None]:
    """Project the explicit Matplotlib embedding payload by position."""
    value = adata.obsm[key]
    coordinates = (
        value.iloc[:, :2].to_numpy(copy=False)
        if isinstance(value, pd.DataFrame)
        else np.asarray(value)[:, :2]
    )
    refs = [obsm(key, 0), obsm(key, 1)]
    xcol, ycol = (plain_name(adata, ref) for ref in refs)
    frame = pd.DataFrame(
        {xcol: coordinates[:, 0], ycol: coordinates[:, 1]},
        index=adata.obs_names.copy(),
    )
    if color is None:
        return frame, xcol, ycol, None

    default_kind, default_layer = resolve_source(adata, layer, use_raw)
    token = parse_token(color)
    cname = plain_name(adata, token)
    if isinstance(token, ObsmRef):
        color_key = embedding_key(adata, token.basis)
        color_value = adata.obsm[color_key]
        if token.index < 0 or token.index >= color_value.shape[1]:
            raise IndexError(
                f"Embedding {color_key!r} has {color_value.shape[1]} coordinates; "
                f"requested index {token.index}."
            )
        vector = (
            color_value.iloc[:, token.index].array
            if isinstance(color_value, pd.DataFrame)
            else np.asarray(color_value)[:, token.index]
        )
    elif isinstance(token, Ref) and token.source == "obs":
        if token.name not in adata.obs.columns:
            raise KeyError(f"obs('{token.name}') not found in adata.")
        vector = adata.obs[token.name].array
    else:
        name = token.name if isinstance(token, Ref) else str(token)
        if not isinstance(token, Ref) and name in adata.obs.columns:
            if name in source_var_names(adata, default_kind):
                warnings.warn(
                    f"'{name}' is both an obs column and a gene; using obs. "
                    f"Use 'gene:{name}' to plot expression.",
                    stacklevel=3,
                )
            vector = adata.obs[name].array
        else:
            kind, selected_layer = default_kind, default_layer
            if isinstance(token, Ref):
                if token.layer is not None:
                    kind, selected_layer = resolve_source(adata, token.layer, token.use_raw)
                elif token.use_raw is not None:
                    kind, selected_layer = resolve_source(adata, None, token.use_raw)
            if name not in source_var_names(adata, kind):
                raise KeyError(f"Could not resolve color={color!r} from obs, genes or obsm.")
            projected, _ = expression_matrix(
                adata,
                [name],
                kind=kind,
                layer=selected_layer,
            )
            vector = _native_vector(projected)
    frame[cname] = vector
    return frame, xcol, ycol, cname


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


def _embedding_style(show_axes: bool | None, equal_aspect: bool | None) -> list:
    """Resolve legacy/publication embedding axes without changing prepared data."""
    publication = _active_style() is not None
    family = "standard" if show_axes is True else "embedding"
    components = [_family_theme(family)]
    if show_axes is None and not publication:
        components.append(_embedding_axes())
    elif show_axes is False:
        components.append(
            theme(
                axis_title=element_blank(),
                axis_text=element_blank(),
                axis_ticks=element_blank(),
                axis_line=element_blank(),
                panel_border=element_blank(),
            )
        )
    if equal_aspect is True or (equal_aspect is None and publication):
        components.append(coord_equal())
    return components


def _centroid_labels(
    df: pd.DataFrame, cname: str, xcol: str, ycol: str, split_by: str | None = None
) -> pd.DataFrame:
    """Median position of each category, for placing a cluster label at its centre.

    When ``split_by`` is given the centroids are computed *within* each facet, so a
    label lands where its category actually sits in that panel rather than at the
    pooled median (which would be broadcast to every facet).
    """
    same = split_by == cname  # the facet column *is* the label column
    keys = [cname] if split_by is None or same else [split_by, cname]
    cents = df.groupby(keys, observed=True)[[xcol, ycol]].median().reset_index()
    cents = cents.rename(columns={cname: "label"})
    if same:
        cents[cname] = cents["label"]
    return cents


_NA_LABEL = "NA"  # scanpy's label for a missing category
_LEGEND_LOCS = {"right margin", "right", "on data", "none"}


def _components_indices(components, width: int, basis: str) -> tuple[int, int]:
    """Zero-based (x, y) coordinate indices from scanpy-style 1-based ``components``."""
    if components is None:
        return 0, 1
    if isinstance(components, str):
        components = components.split(",")
    try:
        pair = tuple(int(c) for c in components)
    except (TypeError, ValueError):
        pair = ()
    if len(pair) != 2 or pair[0] == pair[1] or min(pair) < 1 or max(pair) > width:
        raise ValueError(
            f"components must be two different 1-based coordinates of '{basis}' "
            f"(1..{width}), e.g. (1, 3) or '1,3'; got {components!r}."
        )
    return pair[0] - 1, pair[1] - 1


def _limit(value, values: np.ndarray, name: str):
    """A colour limit: a number, or ``'p<percentile>'`` of the finite values."""
    if value is None:
        return None
    if isinstance(value, str) and value.startswith("p"):
        try:
            return float(np.percentile(values, float(value[1:]))) if len(values) else None
        except ValueError:
            pass
    elif isinstance(value, (int, float, np.integer, np.floating)) and np.isfinite(value):
        return float(value)
    raise ValueError(f"{name} must be a finite number or 'p<percentile>' (e.g. 'p99').")


def _embedding_color_scale(values: pd.Series, *, vmin, vmax, vcenter, color_map, na_color):
    """Continuous colour scale for explicit limits, centring, colour map or NA colour.

    Values outside ``[vmin, vmax]`` are clipped to the end colours, as in scanpy.
    """
    from matplotlib import colormaps
    from matplotlib.colors import TwoSlopeNorm
    from mizani.bounds import squish
    from plotnine import scale_color_gradientn

    finite = values.to_numpy(dtype=float)
    finite = finite[np.isfinite(finite)]
    lo = _limit(vmin, finite, "vmin")
    hi = _limit(vmax, finite, "vmax")
    centre = _limit(vcenter, finite, "vcenter")
    lo = float(finite.min()) if lo is None and len(finite) else lo
    hi = float(finite.max()) if hi is None and len(finite) else hi
    if lo is None or hi is None or lo >= hi:
        raise ValueError(f"vmin ({lo}) must be below vmax ({hi}).")
    cmap = colormaps[color_map or ("RdBu_r" if centre is not None else "viridis")]
    positions = np.linspace(0, 1, 256)
    if centre is None:
        stops = positions
    else:
        if not lo < centre < hi:
            raise ValueError(f"vcenter ({centre}) must lie strictly between vmin and vmax.")
        stops = (TwoSlopeNorm(centre, lo, hi).inverse(positions) - lo) / (hi - lo)
    return scale_color_gradientn(
        colors=[tuple(c) for c in cmap(positions)],
        values=stops,
        limits=(lo, hi),
        oob=squish,
        na_value=na_color,
    )


def _categorical_levels_of(series: pd.Series) -> list:
    if isinstance(series.dtype, pd.CategoricalDtype):
        return list(series.cat.categories)
    return list(pd.Categorical(series).categories)


def _na_labels(breaks) -> list:
    """Legend text for a categorical scale: a missing value reads ``NA``."""
    return [_NA_LABEL if pd.isna(b) or b == "nan" else b for b in breaks]


def _categorical_scale(adata, df: pd.DataFrame, cname: str, *, palette, groups, na_color):
    """Categorical colour scale and the row order that draws missing cells first.

    Missing values stay missing in ``df`` (drawn in ``na_color``, legend ``NA``).
    ``groups`` sets every other category to missing, so only the selected groups
    are coloured and listed. Returns ``(scale, order)``; ``df[cname]`` is
    rewritten in place for ``groups``.
    """
    from mizani.palettes import hue_pal

    from ._palette import _publication_obs_mapping

    series = df[cname]
    levels = _categorical_levels_of(series)
    if groups is not None:
        unknown = [g for g in groups if g not in levels]
        if unknown:
            raise ValueError(f"groups {unknown} are not categories of {cname!r}: {levels}")
        series = df[cname] = series.where(series.isin(groups))
    # missing / de-selected cells first, so coloured cells are never occluded
    order = np.argsort(~series.isna().to_numpy(), kind="stable")
    kwargs = {"na_value": na_color, "labels": _na_labels}
    if groups is not None:
        kwargs["breaks"] = [c for c in levels if c in groups]
    if palette is None and groups is None:
        return scale_color_obs(adata, cname, **kwargs), order
    if isinstance(palette, dict):
        lacking = [c for c in levels if c not in palette]
        if lacking:
            raise ValueError(f"palette has no colour for {lacking}.")
        colours = {c: palette[c] for c in levels}
    elif palette is not None:
        palette = list(palette)
        if len(palette) < len(levels):
            raise ValueError(f"palette has {len(palette)} colours for {len(levels)} categories.")
        colours = dict(zip(levels, palette))
    else:
        mapping = obs_colors(adata, cname) or _publication_obs_mapping(adata, cname)
        colours = (
            {c: mapping[c] for c in levels}
            if mapping is not None
            else dict(zip(levels, hue_pal()(len(levels))))
        )
    return scale_color_manual(values=colours, **kwargs), order


_DENSITY_WARN_CELLS = 20_000


def _kde_params() -> dict:
    """Seed statsmodels' generator so its entropy FutureWarning does not fire.

    The generator only matters for cross-validated bandwidths, which are not
    used; the default ``normal_reference`` bandwidth is deterministic.
    """
    import inspect

    from statsmodels.nonparametric.kernel_density import KDEMultivariate

    if "rng" in inspect.signature(KDEMultivariate.__init__).parameters:
        return {"rng": np.random.default_rng(0)}
    return {}


def _binned_density(xy: np.ndarray) -> np.ndarray:
    """Approximate Gaussian point density: 2-D histogram, Gaussian filter, sample back.

    The bandwidth is statsmodels' ``normal_reference`` rule, which is what the
    exact ``stat_pointdensity`` KDE uses, so the two agree up to binning error.
    Densities are in true density units (integrate to one); each cell reads the
    bilinear interpolation of the smoothed grid at its own position.
    """
    from scipy.ndimage import gaussian_filter, map_coordinates

    n = len(xy)
    spread = xy.std(axis=0)
    if n < 3 or float(spread.min()) == 0.0:
        return np.zeros(n)
    bw = 1.06 * spread * n ** (-1.0 / 6.0)
    # at least ~8 bins per bandwidth so the filter resolves the kernel
    bins = [int(np.clip(np.ptp(xy[:, i]) / bw[i] * 8, 64, 1024)) for i in (0, 1)]
    hist, xe, ye = np.histogram2d(xy[:, 0], xy[:, 1], bins=bins)
    dx, dy = xe[1] - xe[0], ye[1] - ye[0]
    smooth = gaussian_filter(hist, sigma=(bw[0] / dx, bw[1] / dy), mode="constant")
    smooth /= n * dx * dy
    at = [(xy[:, 0] - xe[0]) / dx - 0.5, (xy[:, 1] - ye[0]) / dy - 0.5]
    return map_coordinates(smooth, at, order=1, mode="nearest")


def _point_density(df: pd.DataFrame, xcol: str, ycol: str, split_by: str | None) -> np.ndarray:
    """Binned density per facet group (or all rows), written back by position."""
    xy = df[[xcol, ycol]].to_numpy(dtype=float)
    if split_by is None:
        return _binned_density(xy)
    codes = pd.Categorical(df[split_by]).codes
    density = np.zeros(len(df))
    for code in np.unique(codes):
        rows = np.flatnonzero(codes == code)
        density[rows] = _binned_density(xy[rows])
    return density


def _embedding_panels(adata, basis, colors, titles, ncol, kwargs):
    """One ``plot_embedding`` per colour, composed into a grid (shared option set)."""
    from .layout import compose

    panels = [
        plot_embedding(adata, basis, c, title=t, **kwargs)
        for c, t in zip(colors, titles, strict=True)
    ]
    return compose(panels, ncol=ncol, tag_levels=None)


def plot_embedding(
    adata,
    basis: str = "umap",
    color: str | Sequence[str] | None = None,
    *,
    split_by: str | None = None,
    layer: str | None = None,
    use_raw: bool | None = None,
    components: str | Sequence[int] | None = None,
    size: float = 1.5,
    alpha: float = 0.9,
    pointdensity: bool | None = None,
    method: str = "exact",
    label: bool = False,
    label_size: float = 9,
    low: str = "#d9d9d9",
    high: str = "#2166ac",
    color_map: str | None = None,
    palette: dict | Sequence[str] | None = None,
    vmin: float | str | None = None,
    vmax: float | str | None = None,
    vcenter: float | str | None = None,
    groups: str | Sequence[str] | None = None,
    sort_order: bool = True,
    na_color: str | None = None,
    title: str | None = None,
    legend_loc: str = "right margin",
    ncol: int | None = None,
    downsample: int | None = None,
    random_state: int | None = 0,
    backend: str = "plotnine",
    rasterized: bool = False,
    show_axes: bool | None = None,
    equal_aspect: bool | None = None,
):
    """Scatter over an embedding (UMAP/t-SNE/PCA), optionally coloured and split.

    * ``color=None`` -> density-coloured scatter (``geom_pointdensity``), which
      reads well for dense embeddings.
    * categorical ``color`` (e.g. an obs cluster column) -> discrete colours
      (reusing scanpy's stored palette when present). Missing values are drawn
      in ``na_color`` and labelled ``NA``, as in scanpy.
    * numeric ``color`` (a gene or continuous obs column) -> a gradient.
    * a list of colours -> one panel each, composed into a grid (``ncol``).
    * ``split_by`` -> facet the scatter over an obs column (Seurat ``split.by``).
    * ``label=True`` -> for a categorical ``color``, print each category at its
      centroid using repelled (non-overlapping) text, like scplotter's
      ``CellDimPlot`` / Seurat ``label=TRUE``.
    * ``downsample=N`` -> randomly keep at most ``N`` cells before drawing, for a
      much lighter scatter on large data (density reads the same; exact points differ).
      ``random_state`` controls which cells are retained.

    Parameters
    ----------
    adata
        Annotated data matrix.
    basis : str, default="umap"
        Embedding basis in ``adata.obsm``.
    color : str, selector or sequence, optional
        Observation column, gene, or explicit :func:`ggann.obs`/:func:`ggann.gene`.
        A sequence draws one panel per entry and returns a composition.
    split_by : str, optional
        Observation column used for facets.
    layer, use_raw : optional
        Expression source for gene colour; these options are mutually exclusive.
    components : str or pair of int, optional
        Two different 1-based coordinates to plot, e.g. ``(1, 3)`` or ``"1,3"``
        (scanpy's ``components``). Default ``(1, 2)``.
    size, alpha : float
        Point size and opacity.
    pointdensity : bool, optional
        Use density colouring when ``color`` is absent.
    method : {"exact", "binned"}, default="exact"
        Point-density estimator when ``color`` is absent. ``"exact"`` is the
        Gaussian KDE at every cell (cost grows quadratically; a warning is
        raised above 20,000 cells). ``"binned"`` is an opt-in approximation
        (2-D histogram, Gaussian filter, sampled back at each cell) that is much
        faster; its legend reads ``density (binned)``.
    label : bool, default=False
        Label categorical colour centroids.
    label_size : float
        Centroid-label size.
    low, high : str
        End colours for continuous values.
    color_map : str, optional
        Matplotlib colormap for continuous values (scanpy's ``color_map``);
        overrides ``low``/``high``.
    palette : dict or sequence of str, optional
        Category colours: a ``{category: colour}`` dict or colours in category
        order. Overrides stored ``adata.uns`` colours.
    vmin, vmax, vcenter : float or "p<percentile>", optional
        Continuous colour limits and centre; out-of-range values are clipped to
        the end colours. ``"p99"`` is the 99th percentile of the finite values.
    groups : str or sequence of str, optional
        Categories to colour; all other cells are drawn first in ``na_color``.
    sort_order : bool, default=True
        Draw low values first so high values are not occluded (continuous colour).
    na_color : str, optional
        Colour of missing or de-selected cells. Default: the publication
        missing colour, else light grey.
    title : str, optional
        Plot title. Panels of a colour list default to the colour name.
    legend_loc : {"right margin", "right", "on data", "none"}
        ``"on data"`` is ``label=True`` without a legend; ``"none"`` hides it.
    ncol : int, optional
        Grid columns when ``color`` is a list.
    downsample : int, optional
        Maximum total observations to draw.
    random_state : int, optional
        Reproducible downsampling seed; ``None`` is non-deterministic.
    backend : {"plotnine", "matplotlib"}, default="plotnine"
        Rendering path. The explicit Matplotlib path accelerates an unsplit
        point scatter while returning a composable ggplot subclass. It does not
        support the display options above (``components`` to ``legend_loc``).
    rasterized : bool, default=False
        Rasterize only the point layer. Text, axes, guides, and annotations
        remain editable vectors in SVG and PDF output.
    show_axes : bool, optional
        Show embedding axes when true or hide them when false. The default
        preserves legacy axes outside a publication context and uses the
        publication embedding treatment inside one.
    equal_aspect : bool, optional
        Force equal x/y data units when true. The default enables equal aspect
        only in publication mode.

    Returns
    -------
    plotnine.ggplot
        Composable embedding plot (a composition when ``color`` is a list). The
        Matplotlib backend returns a :class:`ggann.MatplotlibGGPlot` subclass.

    Raises
    ------
    KeyError
        If the embedding or an explicit colour/split source is missing.
    ValueError
        If the embedding has fewer than two coordinates, an option is invalid,
        or the Matplotlib backend is asked for an unsupported density, facet,
        centroid-label or display option.

    Notes
    -----
    Only the requested embedding coordinates and colour metadata are projected.
    Downsampling occurs before extraction and never mutates ``adata``. Cells with
    a missing ``split_by`` are dropped with a warning; cells with a missing
    ``color`` stay visible. Not supported: scanpy's ``projection="3d"``,
    ``edges``/``arrows``/``neighbors_key``, ``add_outline`` and ``frameon``
    (use ``show_axes``); add them with ordinary plotnine layers.

    Examples
    --------
    >>> p = plot_embedding(adata, "umap", color="cell_type")
    >>> p = plot_embedding(adata, "umap", color=["cell_type", "CD3D"], vmax="p99")
    """
    _validate_backend(backend)
    if method not in {"exact", "binned"}:
        raise ValueError("method must be 'exact' or 'binned'.")
    if legend_loc not in _LEGEND_LOCS:
        raise ValueError(f"legend_loc must be one of {sorted(_LEGEND_LOCS)}, got {legend_loc!r}.")
    if isinstance(groups, str):
        groups = [groups]
    adata = _drop_missing_groups(adata, split_by)
    adata = _downsample_cells(adata, None, downsample, random_state=random_state)

    if color is not None and not isinstance(color, (str, Ref, ObsmRef)):
        colors = list(color)
        if not colors:
            raise ValueError("color is an empty sequence.")
        if title is not None and not isinstance(title, str):
            titles = list(title)
        else:
            titles = [title if title is not None else plain_name(adata, c) for c in colors]
        if len(titles) != len(colors):
            raise ValueError(f"{len(titles)} titles for {len(colors)} colours.")
        shared = dict(
            split_by=split_by, layer=layer, use_raw=use_raw, components=components, size=size,
            alpha=alpha, pointdensity=pointdensity, method=method, label=label,
            label_size=label_size, low=low, high=high, color_map=color_map, palette=palette,
            vmin=vmin, vmax=vmax, vcenter=vcenter, groups=groups, sort_order=sort_order,
            na_color=na_color, legend_loc=legend_loc, backend=backend, rasterized=rasterized,
            show_axes=show_axes, equal_aspect=equal_aspect,
        )  # fmt: skip
        return _embedding_panels(adata, basis, colors, titles, ncol, shared)

    if backend == "matplotlib":
        display_options = dict(
            components=components, color_map=color_map, palette=palette, vmin=vmin, vmax=vmax,
            vcenter=vcenter, groups=groups, na_color=na_color, title=title,
        )  # fmt: skip
        used = [k for k, v in display_options.items() if v is not None]
        used += ["sort_order"] * (not sort_order) + ["legend_loc"] * (legend_loc != "right margin")
        used += ["method"] * (method != "exact")
        if used:
            raise ValueError(
                f"backend='matplotlib' does not support {used}; use backend='plotnine'."
            )
        if split_by is not None or label:
            raise ValueError(
                "backend='matplotlib' currently supports unsplit embeddings without "
                "centroid labels; use backend='plotnine' for facets or label=True."
            )
    label_requested = label
    if legend_loc == "on data":
        label = True
    key = embedding_key(adata, basis)
    width = adata.obsm[key].shape[1]
    if width < 2:
        raise ValueError(
            f"Embedding '{basis}' has only {width} dimension(s); "
            "plot_embedding requires at least 2."
        )
    xi, yi = _components_indices(components, width, basis)
    if backend == "matplotlib":
        df, xcol, ycol, native_cname = _native_embedding_frame(
            adata,
            key,
            color,
            layer=layer,
            use_raw=use_raw,
        )
    else:
        coordinate_refs = [obsm(key, xi), obsm(key, yi)]
        requested = [*coordinate_refs]
        if color is not None:
            requested.append(color)
        if split_by is not None:
            requested.append(split_by)
        df = resolve_frame(adata, requested, layer=layer, use_raw=use_raw, strict=True)
        xcol, ycol = (plain_name(adata, ref) for ref in coordinate_refs)
        native_cname = None
    facet = pe.facet_wrap("~" + split_by) if split_by is not None else None
    style = _active_style()
    na_colour = na_color or (style.missing_color if style is not None else "lightgray")
    extras = []
    if title is not None:
        extras.append(labs(title=title))
    if legend_loc in {"none", "on data"}:
        extras.append(theme(legend_position="none"))

    def _with_facet(components):
        return [*components, *extras, *([facet] if facet is not None else [])]

    def _solid_point():
        # plotnine's default 0.5-point same-colour outline duplicates colour
        # conversion and artist work. A zero-width outline with a 0.5 size
        # offset preserves the rendered diameter and represented observations.
        return geom_point(
            size=size + 0.5,
            alpha=alpha,
            stroke=0,
            raster=rasterized,
        )

    def _finish(
        plot,
        *,
        value: str | None = None,
        categorical: bool = False,
        signed: bool = False,
    ):
        if backend == "plotnine":
            return plot
        categories: tuple = ()
        palette: tuple[str, ...] = ()
        if categorical and value is not None:
            stored = obs_colors(adata, value) if value in adata.obs.columns else None
            categories, palette = categorical_palette(df[value], stored)
        return promote_matplotlib_plot(
            plot,
            kind="embedding",
            data=df,
            x=xcol,
            y=ycol,
            value=value,
            categorical=categorical,
            categories=categories,
            palette=palette,
            point_size=size,
            alpha=alpha,
            low=low,
            high=high,
            cmap="Reds",
            publication_style=_active_style(),
            rasterized=rasterized,
            show_axes=show_axes,
            equal_aspect=equal_aspect,
            signed=signed,
        )

    if label_requested and color is None:
        warnings.warn(
            "plot_embedding: label=True is ignored when color is None "
            "(centroid labels need a categorical color).",
            stacklevel=2,
        )

    if color is None:
        if pointdensity is None:
            pointdensity = True
        if pointdensity:
            if backend == "matplotlib":
                raise ValueError(
                    "backend='matplotlib' does not implement point-density "
                    "estimation; pass pointdensity=False or use backend='plotnine'."
                )
            if method == "binned":
                df = df.copy()
                df["density"] = _point_density(df, xcol, ycol, split_by)
                return ggplot(df, aes(xcol, ycol, color="density")) + _with_facet(
                    [
                        geom_point(size=size, alpha=alpha, raster=rasterized),
                        labs(color="density (binned)"),
                        *_embedding_style(show_axes, equal_aspect),
                    ]
                )
            if len(df) > _DENSITY_WARN_CELLS:
                warnings.warn(
                    f"exact point density over {len(df):,} cells is slow (it grows "
                    "quadratically); pass method='binned' for a fast approximation, "
                    "pointdensity=False, or downsample=.",
                    UserWarning,
                    stacklevel=2,
                )
            return ggplot(df, aes(xcol, ycol)) + _with_facet(
                [
                    pe.geom_pointdensity(
                        size=size,
                        alpha=alpha,
                        raster=rasterized,
                        kde_params=_kde_params(),
                    ),
                    labs(color="density"),
                    *_embedding_style(show_axes, equal_aspect),
                ]
            )
        return _finish(
            ggplot(df, aes(xcol, ycol))
            + _with_facet([_solid_point(), *_embedding_style(show_axes, equal_aspect)])
        )

    # `color` may be a bare name, a prefix string ("gene:CD3D@logcounts") or an accessor
    cname = native_cname if native_cname is not None else plain_name(adata, color)
    if cname not in df.columns:
        raise KeyError(f"Could not resolve color={color!r} from obs, genes or obsm.")

    if _is_numeric(df[cname]):
        if label_requested:
            warnings.warn(
                f"plot_embedding: label=True is ignored for the numeric color {color!r} "
                "(centroid labels need a categorical color).",
                stacklevel=2,
            )
        # Draw low-expression cells first so high-expression cells are not occluded
        # (mirrors scanpy's sc.pl.embedding ordering); missing values go first.
        if sort_order:
            df = df.sort_values(cname, na_position="first")
        finite = df[cname].to_numpy(dtype=float)
        finite = finite[np.isfinite(finite)]
        signed = bool(len(finite) and finite.min() < 0 < finite.max())
        custom = any(v is not None for v in (vmin, vmax, vcenter, color_map, na_color))
        if custom and backend == "plotnine":
            continuous_scale = _embedding_color_scale(
                df[cname],
                vmin=vmin,
                vmax=vmax,
                vcenter=vcenter,
                color_map=color_map,
                na_color=na_colour,
            )
        else:
            continuous_scale = (
                _continuous_scale("color", df[cname], "Reds", signed=signed)
                if style is not None
                else scale_color_gradient(low=low, high=high)
            )
        return _finish(
            ggplot(df, aes(xcol, ycol, color=cname))
            + _with_facet(
                [
                    _solid_point(),
                    continuous_scale,
                    *_embedding_style(show_axes, equal_aspect),
                ]
            ),
            value=cname,
            signed=signed,
        )
    colour_scale = scale_color_obs(adata, cname)
    if backend == "plotnine":
        df = df.copy()
        colour_scale, order = _categorical_scale(
            adata, df, cname, palette=palette, groups=groups, na_color=na_colour
        )
        df = df.iloc[order]
    panel = [
        _solid_point(),
        colour_scale,
        # enlarge the legend swatches so categories stay readable (scplotter does this)
        guides(
            color=guide_legend(
                override_aes={"size": 4},
                ncol=(2 if style is not None and int(df[cname].nunique()) > 8 else None),
                byrow=style is not None,
            )
        ),
        *_embedding_style(show_axes, equal_aspect),
    ]
    if label:
        cents = _centroid_labels(df[df[cname].notna()], cname, xcol, ycol, split_by=split_by)
        # white-backed repelled labels at centroids, like scplotter's label_bg="white"
        panel.append(
            pe.geom_label_repel(
                aes(xcol, ycol, label="label"),
                data=cents,
                size=label_size * 0.85,
                fill="white",
                color="black",
                inherit_aes=False,
            )
        )
    return _finish(
        ggplot(df, aes(xcol, ycol, color=cname)) + _with_facet(panel),
        value=cname,
        categorical=True,
    )


def plot_features(
    adata,
    features: Sequence[str],
    basis: str = "umap",
    *,
    layer: str | None = None,
    use_raw: bool | None = None,
    components: str | Sequence[int] | None = None,
    ncol: int | None = None,
    size: float = 1.2,
    alpha: float = 0.9,
    cmap: str = "magma",
    vmin: float | str | None = None,
    vmax: float | str | None = None,
    vcenter: float | str | None = None,
    sort_order: bool = True,
    na_color: str | None = None,
    color_scale: str = "auto",
    downsample: int | None = None,
    random_state: int | None = 0,
    rasterized: bool = False,
    show_axes: bool | None = None,
    equal_aspect: bool | None = None,
):
    """Multi-gene embedding grid: one panel per feature.

    ``color_scale="shared"`` draws one faceted plot whose panels share a single
    expression colour scale, best for comparing magnitudes *across* genes of one
    unit. ``"free"`` draws one :func:`plot_embedding` per feature, each with its
    own colour bar and original values, composed into a grid. ``"auto"`` (the
    default) is ``"free"`` as soon as a feature is an obs column (counts,
    scores and expression are different units) and ``"shared"`` for genes only.
    ``downsample=N`` caps cells before drawing; ``random_state`` controls them.

    Parameters
    ----------
    adata
        Annotated data matrix.
    features : sequence of str
        Genes or numeric observation columns.
    basis : str, default="umap"
        Embedding basis.
    layer, use_raw : optional
        Mutually exclusive expression-source selection.
    components : str or pair of int, optional
        Two different 1-based coordinates to plot, e.g. ``(1, 3)``.
    ncol : int, optional
        Panel columns.
    size, alpha : float
        Point size and opacity.
    cmap : str
        Matplotlib colormap name.
    vmin, vmax, vcenter : float or "p<percentile>", optional
        Colour limits and centre; out-of-range values clip to the end colours.
        In the shared mode they apply to all panels at once.
    sort_order : bool, default=True
        Draw low values first so high values are not occluded.
    na_color : str, optional
        Colour of missing values.
    color_scale : {"auto", "shared", "free"}
        Shared or per-feature colour scales (see above).
    downsample : int, optional
        Maximum total observations to draw.
    random_state : int, optional
        Downsampling seed.
    rasterized : bool, default=False
        Rasterize only the embedding points in vector output.
    show_axes : bool, optional
        Explicitly show or hide embedding axes; ``None`` selects the active
        legacy or publication default.
    equal_aspect : bool, optional
        Force equal x/y data units; ``None`` enables it in publication mode.

    Returns
    -------
    plotnine.ggplot
        Composable faceted feature plot; a composition for ``color_scale="free"``.

    Raises
    ------
    KeyError
        If the embedding or an explicit source is missing.
    ValueError
        If no numeric feature resolves, the embedding is one-dimensional, or
        downsampling/source selection is invalid.

    Notes
    -----
    annplyr projects requested fields only; ``adata`` is never mutated. Values
    are never rescaled; per-feature normalization is not offered here.

    Examples
    --------
    >>> p = plot_features(adata, ["CD3D", "NKG7"], basis="umap")
    >>> p = plot_features(adata, ["CD3D", "n_counts"], color_scale="free")
    """
    features = ordered_unique(features)
    adata = _downsample_cells(adata, None, downsample, random_state=random_state)
    if color_scale not in {"auto", "shared", "free"}:
        raise ValueError("color_scale must be 'auto', 'shared' or 'free'.")
    key = embedding_key(adata, basis)
    width = adata.obsm[key].shape[1]
    if width < 2:
        raise ValueError(f"Embedding '{basis}' has fewer than 2 dimensions.")
    xi, yi = _components_indices(components, width, basis)
    coordinate_refs = [obsm(key, xi), obsm(key, yi)]
    coords = resolve_frame(adata, coordinate_refs)
    xcol, ycol = coords.columns[:2]

    values = resolve_frame(adata, list(features), layer=layer, use_raw=use_raw, strict=True)
    if not features:
        raise ValueError("plot_features needs at least one numeric feature (gene or metric).")
    feats = [plain_name(adata, f) for f in features]
    non_numeric = [f for f in feats if not _is_numeric(values[f])]
    if non_numeric:
        raise ValueError(
            f"plot_features plots numeric features only; {non_numeric} are categorical. "
            "Use plot_embedding(color=...) for those."
        )

    mixed_units = any(f in adata.obs.columns for f in feats)
    if color_scale == "free" or (color_scale == "auto" and mixed_units):
        return _embedding_panels(
            adata,
            basis,
            list(features),
            [str(f) for f in feats],
            ncol,
            dict(
                layer=layer, use_raw=use_raw, components=components, size=size, alpha=alpha,
                color_map=cmap, vmin=vmin, vmax=vmax, vcenter=vcenter, sort_order=sort_order,
                na_color=na_color, rasterized=rasterized, show_axes=show_axes,
                equal_aspect=equal_aspect,
            ),
        )  # fmt: skip

    df = coords.join(values[feats])
    long = df.melt(
        id_vars=[xcol, ycol],
        value_vars=feats,
        var_name="feature",
        value_name="expression",
    )
    if sort_order:
        long = long.sort_values("expression", na_position="first")
    long["feature"] = pd.Categorical(long["feature"], categories=feats, ordered=True)
    if any(v is not None for v in (vmin, vmax, vcenter, na_color)):
        style = _active_style()
        scale = _embedding_color_scale(
            long["expression"],
            vmin=vmin,
            vmax=vmax,
            vcenter=vcenter,
            color_map=cmap,
            na_color=na_color or (style.missing_color if style is not None else "lightgray"),
        )
    else:
        scale = _continuous_scale("color", long["expression"], cmap, signed=False)
    return (
        ggplot(long, aes(xcol, ycol, color="expression"))
        + geom_point(size=size, alpha=alpha, raster=rasterized)
        + pe.facet_wrap("~feature", ncol=ncol)
        + scale
        + _embedding_style(show_axes, equal_aspect)
    )


def _cell_rank(tidy: pd.DataFrame, group_by: str) -> pd.DataFrame:
    """Add a ``cell_rank`` column that orders cells by their ``group_by`` category.

    Cells are sorted by group then numbered 0..N-1, so a per-cell x layout reads
    left-to-right by group. Shared by :func:`plot_heatmap` and
    ``markers.plot_tracksplot``.
    """
    positions = tidy.groupby("feature", observed=True, sort=False).cumcount()
    cell = cast(
        pd.DataFrame,
        tidy.loc[:, ["obs_name", group_by]].drop_duplicates(),
    )
    n_cells = int(positions.max()) + 1 if len(positions) else 0
    if cell["obs_name"].is_unique and len(cell) == n_cells:
        cell = cell.sort_values(group_by, kind="stable")
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


_MAX_MATRIX_GROUPS = 100  # beyond this a dotplot/matrixplot is unreadable
_KEY_SEP = "\x1f"  # joins marker group and gene into one unique axis key


def _marker_genes(genes):
    """Split ``genes`` into ``(unique genes, occurrences)``.

    ``genes`` is a list, or a ``{marker group: [genes]}`` dict as in scanpy. For
    a dict ``occurrences`` lists every ``(group, gene)`` pair, so a gene listed
    in two groups is drawn twice; its values come from one computation.
    ``occurrences`` is ``None`` for a plain list.
    """
    if not isinstance(genes, dict):
        return ordered_unique(genes), None
    occurrences = [
        (str(label), gene) for label, members in genes.items() for gene in ordered_unique(members)
    ]
    if not occurrences:
        raise ValueError("var_names dict contains no genes.")
    return ordered_unique(gene for _, gene in occurrences), occurrences


def _expand_occurrences(frame: pd.DataFrame, occurrences, *, reverse: bool = False):
    """Repeat each gene's rows once per marker-group occurrence.

    Adds ``marker_group`` and ``feature_key`` (a unique, display-only axis key);
    ``reverse`` orders the keys last-first for an axis that draws bottom-up.
    """
    parts = []
    for label, gene in occurrences:
        part = frame[frame["feature"] == gene].copy()
        part["marker_group"] = label
        part["feature_key"] = f"{label}{_KEY_SEP}{gene}"
        parts.append(part)
    out = pd.concat(parts, ignore_index=True)
    keys = [f"{label}{_KEY_SEP}{gene}" for label, gene in occurrences]
    out["feature_key"] = pd.Categorical(
        out["feature_key"], categories=keys[::-1] if reverse else keys
    )
    out["marker_group"] = pd.Categorical(
        out["marker_group"], categories=list(dict.fromkeys(label for label, _ in occurrences))
    )
    return out


def _key_labels(keys) -> list[str]:
    return [str(key).split(_KEY_SEP, 1)[-1] for key in keys]


def _colour_limits(scale, values: pd.Series, vmin, vmax, *, symmetric: bool = False):
    """Set explicit (or zero-centred) limits on a continuous scale; clip outliers.

    Returns the numeric ``(vmin, vmax)`` actually used, or ``(None, None)``.
    """
    if vmin is None and vmax is None and not symmetric:
        return None, None
    from mizani.bounds import squish

    finite = values.to_numpy(dtype=float)
    finite = finite[np.isfinite(finite)]
    low, high = _limit(vmin, finite, "vmin"), _limit(vmax, finite, "vmax")
    if symmetric and low is None and high is None and len(finite):
        high = float(np.abs(finite).max()) or 1.0
        low = -high
    low = float(finite.min()) if low is None and len(finite) else low
    high = float(finite.max()) if high is None and len(finite) else high
    if low is None or high is None or low >= high:
        raise ValueError(f"vmin ({low}) must be below vmax ({high}).")
    scale.limits = (low, high)
    scale.oob = squish
    return low, high


def _matrix_layout(
    group_by, group_levels, genes, occurrences, split_by, swap_axes
) -> tuple[str, str, list]:
    """Axis mapping, scales and facets of a groups-by-features display.

    The first group (or gene when swapped) is drawn at the top, as in scanpy.
    Only the displayed axis is reversed: data, aggregation and category order
    keep their semantic order. ``occurrences`` (from a marker dict) facet the
    feature axis into one panel per marker group.
    """
    feature = "feature_key" if occurrences else "feature"
    xcol, ycol = (group_by, feature) if swap_axes else (feature, group_by)
    parts: list = []
    if not swap_axes:
        parts.append(scale_y_discrete(limits=list(reversed(group_levels))))
    elif not occurrences:
        parts.append(scale_y_discrete(limits=list(reversed(genes))))
    if occurrences:  # keys are pre-reversed by _expand_occurrences when on y
        labels = {"labels": _key_labels}
        parts.append(scale_y_discrete(**labels) if swap_axes else scale_x_discrete(**labels))
        rows = "marker_group" if swap_axes else (split_by or ".")
        cols = (split_by or ".") if swap_axes else "marker_group"
        free = "free_y" if swap_axes else "free_x"
        parts.append(facet_grid(f"{rows} ~ {cols}", scales=free, space=free))
    elif split_by:
        parts.append(pe.facet_wrap(f"~{split_by}"))
    return xcol, ycol, parts


def plot_dotplot(
    adata,
    genes: Sequence[str] | dict[str, Sequence[str]],
    group_by: str,
    *,
    split_by: str | None = None,
    layer: str | None = None,
    use_raw: bool | None = None,
    standard_scale: str | None = None,
    expression_cutoff: float = 0.0,
    cmap: str | None = None,
    vmin: float | str | None = None,
    vmax: float | str | None = None,
    dot_min: float | None = None,
    dot_max: float | None = None,
    swap_axes: bool = False,
    size_range: tuple[float, float] = (0.5, 8.0),
    categories_order: Iterable[str] | None = None,
    backend: str = "plotnine",
    rasterized: bool = False,
    annotate: bool | str = False,
):
    """Marker dotplot: dot *size* = fraction expressing, *colour* = mean expression.

    Defaults to ``adata.raw`` so values match ``sc.pl.dotplot``. ``split_by`` adds a
    facet column so the dotplot is repeated per split level (scplotter ``split_by``).

    Parameters
    ----------
    adata
        Annotated data matrix.
    genes : sequence of str or dict
        Genes to summarize, in display order. A ``{marker group: [genes]}``
        dict (scanpy's ``var_names`` dict) draws one panel per marker group; a
        gene listed in several groups appears once per occurrence, labelled by
        group, with values from a single computation. Needs ``backend="plotnine"``.
    group_by : str
        Primary observation grouping column.
    split_by : str, optional
        Additional grouping and facet column.
    layer, use_raw : optional
        Mutually exclusive expression-source selection.
    standard_scale : {None, "var", "group", "zscore"}
        Optional scaling of summarized means.
    expression_cutoff : float
        Threshold used for the expressing fraction.
    cmap : str, optional
        Matplotlib colormap name. Default ``"Reds"``, or the diverging
        ``"RdBu_r"`` with zero-centred limits for ``standard_scale="zscore"``.
    vmin, vmax : float or "p<percentile>", optional
        Colour limits; values outside are clipped to the end colours.
    dot_min, dot_max : float, optional
        Fractions (0 to 1) mapped to the smallest and largest dot; others clip.
    swap_axes : bool, default=False
        Put genes on y and groups on x. The first group (or gene) is always
        drawn at the top; use this rather than ``coord_flip()``, which would
        reverse the displayed order.
    size_range : tuple of float
        Minimum and maximum dot size.
    categories_order : iterable of str, optional
        Complete order of observed groups.
    backend : {"plotnine", "matplotlib"}, default="plotnine"
        Rendering path. The explicit Matplotlib path supports the unsplit layout
        and returns a composable ggplot subclass.
    rasterized : bool, default=False
        Rasterize the dot artist while keeping labels and guides as vectors.
    annotate : bool or {"auto", "force"}, default=False
        Add mean-expression labels. ``"auto"`` draws a label only when the
        rendered cell is at least 12 points wide and high; ``"force"`` always
        draws labels. ``True`` is an alias for ``"force"``.

    Returns
    -------
    plotnine.ggplot
        Composable dotplot; the explicit backend returns
        :class:`ggann.MatplotlibGGPlot`.

    Raises
    ------
    KeyError
        If a gene, grouping column, or selected layer is missing.
    ValueError
        If source, scaling, category ordering, backend, or the requested direct
        layout is invalid.

    Notes
    -----
    Requested genes are projected before sparse-native aggregation; no implicit
    downsampling or whole-matrix materialization occurs. ``adata`` is unchanged.
    A warning is raised when ``group_by`` has more than 100 levels (a numeric
    column is one level per value). Not supported: scanpy's ``dendrogram``,
    ``var_group_rotation``, ``smallest_dot`` and ``mean_only_expressed``.

    Examples
    --------
    >>> p = plot_dotplot(adata, ["CD3D", "NKG7"], group_by="cell_type")
    >>> p = plot_dotplot(adata, {"T": ["CD3D"], "NK": ["NKG7"]}, "cell_type")
    """
    adata = _drop_missing_groups(adata, group_by, split_by)
    _check_group_cardinality(adata, group_by, _MAX_MATRIX_GROUPS, action="warn")
    _validate_backend(backend)
    if backend == "matplotlib" and (split_by is not None or isinstance(genes, dict)):
        raise ValueError(
            "backend='matplotlib' currently supports unsplit dotplots of a gene list; "
            "use backend='plotnine' for split_by facets or a marker dict."
        )
    genes, occurrences = _marker_genes(genes)
    zscore = standard_scale == "zscore"
    cmap = cmap or ("RdBu_r" if zscore else "Reds")
    extra = [split_by] if split_by else []
    aggregate = aggregate_expression_native if backend == "matplotlib" else aggregate_expression
    agg = aggregate(
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
    if occurrences:
        agg = _expand_occurrences(agg, occurrences, reverse=swap_axes)
    xcol, ycol, layout = _matrix_layout(
        group_by, list(agg[group_by].cat.categories), genes, occurrences, split_by, swap_axes
    )
    color_lab = "scaled\nexpression" if standard_scale else "mean\nexpression"
    size_kwargs = {
        "range": size_range,
        "labels": lambda xs: [f"{x:.0%}" for x in xs],
    }
    if _active_style() is not None:
        size_kwargs.update(limits=(0, 1), breaks=(0, 0.25, 0.5, 0.75, 1.0))
    dot_limits = None
    if dot_min is not None or dot_max is not None:
        from mizani.bounds import squish

        lo = 0.0 if dot_min is None else float(dot_min)
        hi = float(agg["fraction"].max()) if dot_max is None else float(dot_max)
        if not 0 <= lo < hi:
            raise ValueError(f"dot_min ({lo}) must be non-negative and below dot_max ({hi}).")
        dot_limits = (lo, hi)
        size_kwargs.update(limits=dot_limits, oob=squish)
    colour_scale = _continuous_scale("color", agg["mean_expression"], cmap, signed=zscore)
    limits = _colour_limits(
        colour_scale,
        agg["mean_expression"],
        vmin,
        vmax,
        symmetric=zscore and _active_style() is None,
    )
    components = [
        geom_point(
            aes(size="fraction", color="mean_expression"),
            raster=rasterized,
        ),
        colour_scale,
        scale_size(**size_kwargs),
        labs(x="", y="", color=color_lab, size="fraction\nexpressing"),
        _family_theme("matrix"),
        pe.rotate_x_text(45),
        *layout,
    ]
    threshold = annotation_threshold(annotate)
    if threshold is not None:
        style = _active_style()
        components.append(
            geom_contrast_text(
                aes(label="mean_expression"),
                min_cell_pt=threshold,
                format_string="{:.2g}",
                size=style.axis_text_size if style is not None else 7,
            )
        )
    plot = ggplot(agg, aes(xcol, ycol)) + components
    if backend == "plotnine":
        return plot
    return promote_matplotlib_plot(
        plot,
        kind="dotplot",
        data=agg,
        x=xcol,
        y=ycol,
        value="mean_expression",
        fraction="fraction",
        cmap=cmap,
        size_range=size_range,
        value_label=color_lab,
        publication_style=_active_style(),
        rasterized=rasterized,
        signed=zscore,
        annotation_min_cell_pt=threshold,
        vmin=limits[0],
        vmax=limits[1],
        symmetric=zscore,
        dot_limits=dot_limits,
    )


def plot_matrixplot(
    adata,
    genes: Sequence[str] | dict[str, Sequence[str]],
    group_by: str,
    *,
    split_by: str | None = None,
    layer: str | None = None,
    use_raw: bool | None = None,
    standard_scale: str | None = None,
    cmap: str | None = None,
    vmin: float | str | None = None,
    vmax: float | str | None = None,
    swap_axes: bool = False,
    categories_order: Iterable[str] | None = None,
    backend: str = "plotnine",
    rasterized: bool = False,
    annotate: bool | str = False,
):
    """Aggregated mean-expression heatmap (genes x groups) as a plotnine tile plot.

    Like ``sc.pl.matrixplot``, this defaults to ``standard_scale=None`` (raw group
    means). Pass ``standard_scale="var"`` to rescale each gene to ``[0, 1]``, which
    keeps cross-gene patterns legible when magnitudes differ widely. ``split_by``
    adds a facet column so the heatmap is repeated per split level.

    Parameters
    ----------
    adata
        Annotated data matrix.
    genes : sequence of str or dict
        Genes to summarize, in display order. A ``{marker group: [genes]}``
        dict draws one panel per marker group; a gene listed in several groups
        appears once per occurrence, labelled by group, with values from a single
        computation. Needs ``backend="plotnine"``.
    group_by : str
        Primary observation grouping column.
    split_by : str, optional
        Additional grouping and facet column.
    layer, use_raw : optional
        Mutually exclusive expression-source selection.
    standard_scale : {None, "var", "group", "zscore"}
        Optional scaling of group means.
    cmap : str, optional
        Matplotlib colormap name. Default ``"viridis"``, or the diverging
        ``"RdBu_r"`` with zero-centred limits for ``standard_scale="zscore"``.
    vmin, vmax : float or "p<percentile>", optional
        Colour limits; values outside are clipped to the end colours.
    swap_axes : bool, default=False
        Put genes on y and groups on x. The first group (or gene) is always
        drawn at the top; use this rather than ``coord_flip()``.
    categories_order : iterable of str, optional
        Complete order of observed groups.
    backend : {"plotnine", "matplotlib"}, default="plotnine"
        Rendering path. The explicit Matplotlib path supports the unsplit layout
        and returns a composable ggplot subclass.
    rasterized : bool, default=False
        Rasterize the tile layer while keeping labels and guides as vectors.
    annotate : bool or {"auto", "force"}, default=False
        Add mean-expression labels. ``"auto"`` requires each rendered cell to
        be at least 12 points in both dimensions; ``"force"`` always draws.

    Returns
    -------
    plotnine.ggplot
        Composable mean-expression tile plot; the explicit backend returns
        :class:`ggann.MatplotlibGGPlot`.

    Raises
    ------
    KeyError
        If a gene, grouping column, or selected layer is missing.
    ValueError
        If source, scaling, category ordering, backend, or the requested direct
        layout is invalid.

    Notes
    -----
    Requested genes are projected before sparse-native aggregation. ``adata`` is
    not mutated. A warning is raised when ``group_by`` has more than 100 levels.
    A matrixplot has no detection threshold, so ``expression_cutoff`` and the
    ``dot_*`` options belong to :func:`plot_dotplot` only. Not supported:
    scanpy's ``dendrogram`` and ``var_group_rotation``.

    Examples
    --------
    >>> p = plot_matrixplot(adata, ["CD3D", "NKG7"], group_by="cell_type")
    """
    adata = _drop_missing_groups(adata, group_by, split_by)
    _check_group_cardinality(adata, group_by, _MAX_MATRIX_GROUPS, action="warn")
    _validate_backend(backend)
    if backend == "matplotlib" and (split_by is not None or isinstance(genes, dict)):
        raise ValueError(
            "backend='matplotlib' currently supports unsplit matrixplots of a gene list; "
            "use backend='plotnine' for split_by facets or a marker dict."
        )
    genes, occurrences = _marker_genes(genes)
    zscore = standard_scale == "zscore"
    cmap = cmap or ("RdBu_r" if zscore else "viridis")
    extra = [split_by] if split_by else []
    aggregate = aggregate_means_native if backend == "matplotlib" else aggregate_means
    agg = aggregate(
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
    if occurrences:
        agg = _expand_occurrences(agg, occurrences, reverse=swap_axes)
    xcol, ycol, layout = _matrix_layout(
        group_by, list(agg[group_by].cat.categories), genes, occurrences, split_by, swap_axes
    )
    color_lab = "scaled\nexpression" if standard_scale else "mean\nexpression"
    fill_scale = _continuous_scale("fill", agg["mean_expression"], cmap, signed=zscore)
    limits = _colour_limits(
        fill_scale,
        agg["mean_expression"],
        vmin,
        vmax,
        symmetric=zscore and _active_style() is None,
    )
    components = [
        geom_tile(raster=rasterized),
        fill_scale,
        labs(x="", y="", fill=color_lab),
        _family_theme("matrix"),
        pe.rotate_x_text(45),
        *layout,
    ]
    threshold = annotation_threshold(annotate)
    if threshold is not None:
        style = _active_style()
        components.append(
            geom_contrast_text(
                aes(label="mean_expression"),
                min_cell_pt=threshold,
                format_string="{:.2g}",
                size=style.axis_text_size if style is not None else 7,
            )
        )
    plot = ggplot(agg, aes(xcol, ycol, fill="mean_expression")) + components
    if backend == "plotnine":
        return plot
    return promote_matplotlib_plot(
        plot,
        kind="matrixplot",
        data=agg,
        x=xcol,
        y=ycol,
        value="mean_expression",
        cmap=cmap,
        value_label=color_lab,
        publication_style=_active_style(),
        rasterized=rasterized,
        signed=zscore,
        annotation_min_cell_pt=threshold,
        vmin=limits[0],
        vmax=limits[1],
        symmetric=zscore,
    )


def plot_embedding_density(
    adata,
    basis: str = "umap",
    group_by: str | None = None,
    *,
    groups: str | Sequence[str] | None = None,
    background: bool | None = None,
    components: str | Sequence[int] | None = None,
    size: float = 2.0,
    ncol: int | None = None,
    cmap: str = "viridis",
    vmin: float | str | None = None,
    vmax: float | str | None = None,
    downsample: int | None = None,
    random_state: int | None = 0,
    rasterized: bool = False,
    show_axes: bool | None = None,
    equal_aspect: bool | None = None,
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

    Parameters
    ----------
    adata
        Annotated data matrix.
    basis : str, default="umap"
        Embedding basis.
    group_by : str, optional
        Observation column defining density facets. At most 100 levels; a
        numeric column with more raises ``ValueError`` before any density is
        computed.
    groups : str or sequence of str, optional
        Draw only these groups (scanpy's ``group``). Each group's density is
        the same as when all groups are drawn.
    background : bool, optional
        Draw every cell in grey behind each panel, as scanpy does. The grey
        cells are context only: they are not rows of ``plot.data`` and take no
        part in any density. Default: on when ``groups`` is given, else off.
    components : str or pair of int, optional
        Two different 1-based coordinates to plot, e.g. ``(1, 3)``.
    size : float
        Point size.
    ncol : int, optional
        Facet columns.
    cmap : str
        Matplotlib colormap name.
    vmin, vmax : float or "p<percentile>", optional
        Colour limits on the ``[0, 1]`` density; values outside clip.
    downsample : int, optional
        Maximum observations per group.
    random_state : int, optional
        Downsampling seed.
    rasterized : bool, default=False
        Rasterize only density points in vector output.
    show_axes : bool, optional
        Explicitly show or hide embedding axes; ``None`` selects the active
        legacy or publication default.
    equal_aspect : bool, optional
        Force equal x/y data units; ``None`` enables it in publication mode.

    Returns
    -------
    plotnine.ggplot
        Composable density plot.

    Raises
    ------
    KeyError
        If the embedding or grouping column is missing.
    ValueError
        If the embedding is one-dimensional or downsampling is invalid.

    Notes
    -----
    KDE cost grows with retained cells; downsampling is explicit and deterministic
    by default. The input is not mutated. Densities are computed here, never read
    from ``sc.tl.embedding_density`` results (scanpy's ``key`` is not supported).

    Examples
    --------
    >>> p = plot_embedding_density(adata, "umap", group_by="cell_type")
    """
    adata = _drop_missing_groups(adata, group_by)
    _check_group_cardinality(adata, group_by, _MAX_DENSITY_GROUPS, action="raise")
    if isinstance(groups, str):
        groups = [groups]
    if groups is not None and group_by is None:
        raise ValueError("groups needs group_by.")
    adata = _downsample_cells(adata, group_by, downsample, random_state=random_state)
    key = embedding_key(adata, basis)
    width = adata.obsm[key].shape[1]
    if width < 2:
        raise ValueError(f"Embedding '{basis}' has fewer than 2 dimensions.")
    xi, yi = _components_indices(components, width, basis)
    coords = resolve_frame(adata, [obsm(key, xi), obsm(key, yi)])
    xcol, ycol = coords.columns[:2]
    xy = coords[[xcol, ycol]].to_numpy(dtype=float)

    if group_by is None:
        df = coords.copy()
        df["density"] = _kde_density(xy.T)
        plot = ggplot(df, aes(xcol, ycol, color="density"))
    else:
        # Everything is positional: a group's density is written back to the
        # rows of its own cells, so interleaved groups and duplicate obs names
        # cannot misalign (and a single group is no special case).
        labels = pd.Categorical(adata.obs[group_by])
        cats = list(labels.categories)
        levels = [str(c) for c in cats]
        if groups is not None:
            unknown = [g for g in groups if str(g) not in levels]
            if unknown:
                raise ValueError(f"groups {unknown} are not levels of {group_by!r}: {levels}")
            wanted = {str(g) for g in groups}
        else:
            wanted = set(levels)
        density = np.zeros(len(coords))
        selected = np.zeros(len(coords), dtype=bool)
        for code, level in enumerate(levels):
            if level not in wanted:
                continue
            rows = np.flatnonzero(labels.codes == code)
            density[rows] = _kde_density(xy[rows].T)
            selected[rows] = True
        df = coords.copy()
        df[group_by] = adata.obs[group_by].array
        df["density"] = density
        df = df[selected].copy()
        order = [c for c in cats if str(c) in wanted]
        df = _order_groups(df, group_by, order)
        plot = ggplot(df, aes(xcol, ycol, color="density"))
        if background is None:
            background = groups is not None
        if background:
            # grey context in every panel; layer data, not observations of plot.data
            plot = plot + geom_point(
                aes(xcol, ycol),
                data=coords[[xcol, ycol]],
                color="#d9d9d9",
                size=size,
                raster=rasterized,
                inherit_aes=False,
            )
    plot = plot + geom_point(size=size, raster=rasterized)
    if group_by is not None:
        plot = plot + pe.facet_wrap("~" + group_by, ncol=ncol)
    scale = _continuous_scale("color", df["density"], cmap, signed=False)
    _colour_limits(scale, df["density"], vmin, vmax)
    return plot + scale + _embedding_style(show_axes, equal_aspect)


def plot_heatmap(
    adata,
    genes: Sequence[str] | dict[str, Sequence[str]],
    group_by: str,
    *,
    layer: str | None = None,
    use_raw: bool | None = None,
    cmap: str | None = None,
    vmin: float | str | None = None,
    vmax: float | str | None = None,
    swap_axes: bool = False,
    standard_scale: str | None = None,
    categories_order: Sequence[str] | None = None,
    downsample: int | None = None,
    random_state: int | None = 0,
    rasterized: bool = False,
    annotate: bool | str = False,
):
    """Per-**cell** expression heatmap, cells grouped along x (``sc.pl.heatmap``).

    One column per cell (blocked and labelled by ``group_by``), one row per gene,
    tile-coloured by expression -- the per-cell counterpart to the aggregated
    :func:`plot_matrixplot`. ``standard_scale='var'`` z-min-max-scales each gene to
    ``[0, 1]`` across cells so low- and high-range genes stay comparable.
    ``downsample=N`` caps cells per group first; ``random_state`` controls which
    cells are retained.

    Parameters
    ----------
    adata
        Annotated data matrix.
    genes : sequence of str or dict
        Genes shown as rows. A ``{marker group: [genes]}`` dict adds one row
        panel per marker group; a gene listed in several groups appears once per
        occurrence, labelled by group, with values from a single computation.
    group_by : str
        Observation column defining cell blocks.
    layer, use_raw : optional
        Mutually exclusive expression-source selection.
    cmap : str, optional
        Matplotlib colormap name. Default ``"viridis"``, or the diverging
        ``"RdBu_r"`` with zero-centred limits for ``standard_scale="zscore"``.
    vmin, vmax : float or "p<percentile>", optional
        Colour limits; values outside are clipped to the end colours.
    swap_axes : bool, default=False
        Draw genes as columns and cells as rows (groups stacked top to bottom).
    standard_scale : {None, "var", "group", "zscore"}
        Optional per-gene or per-cell scaling.
    categories_order : sequence of str, optional
        Complete order of observed groups.
    downsample : int, optional
        Maximum observations retained per group.
    random_state : int, optional
        Downsampling seed.
    rasterized : bool, default=False
        Rasterize the expression tiles while preserving vector text and axes.
    annotate : bool or {"auto", "force"}, default=False
        Add expression labels. ``"auto"`` draws only in cells at least 12
        points wide and high; ``"force"`` always draws.

    Returns
    -------
    plotnine.ggplot
        Composable per-cell heatmap.

    Raises
    ------
    KeyError
        If a gene, grouping column, or layer is missing.
    ValueError
        If source, scaling, ordering, or downsampling is invalid.

    Notes
    -----
    The prepared long table has one row per retained cell and gene; use explicit
    downsampling for large inputs. ``adata`` is not mutated.

    Examples
    --------
    >>> p = plot_heatmap(adata, ["CD3D", "NKG7"], group_by="cell_type")
    """
    adata = _drop_missing_groups(adata, group_by)
    adata = _downsample_cells(adata, group_by, downsample, random_state=random_state)
    genes, occurrences = _marker_genes(genes)
    zscore = standard_scale == "zscore"
    cmap = cmap or ("RdBu_r" if zscore else "viridis")
    tidy = tidy_expression(adata, genes, group_by, layer=layer, use_raw=use_raw)
    if categories_order is None:
        categories_order = _group_categories(adata, group_by)
    tidy = _order_groups(tidy, group_by, categories_order)

    if standard_scale not in {None, "var", "group", "zscore"}:
        raise ValueError(
            f"standard_scale must be None, 'var', 'group' or 'zscore', got {standard_scale!r}"
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
            rng = (g.transform("max") - lo).replace(0, 1)  # pyright: ignore[reportOperatorIssue]
            tidy["value"] = (tidy["value"] - lo) / rng
        fill_lab = "scaled expr."
    else:
        fill_lab = "expression"

    # order cells by group, then give each a rank so tiles sit side by side
    tidy = _cell_rank(tidy, group_by)
    if occurrences:
        tidy = _expand_occurrences(tidy, occurrences, reverse=not swap_axes)
        feature = "feature_key"
    else:
        feature = "feature"
        # first gene on top (y) or left (x)
        tidy["feature"] = pd.Categorical(
            tidy["feature"],
            categories=genes if swap_axes else list(reversed(genes)),
            ordered=True,
        )

    marker = "marker_group" if occurrences else "."
    fill_scale = _continuous_scale("fill", tidy["value"], cmap, signed=zscore)
    _colour_limits(
        fill_scale,
        tidy["value"],
        vmin,
        vmax,
        symmetric=zscore and _active_style() is None,
    )
    if swap_axes:
        xcol, ycol = feature, "cell_rank"
        layout = [
            facet_grid(f"{group_by} ~ {marker}", scales="free", space="free"),
            scale_y_reverse(),
        ]
        hide = theme(
            axis_text_y=element_blank(),
            axis_ticks_major_y=element_blank(),
            strip_text_y=element_text(angle=0),
            strip_text_x=element_text(angle=0),
        )
    else:
        xcol, ycol = "cell_rank", feature
        layout = [facet_grid(f"{marker} ~ {group_by}", scales="free", space="free")]
        hide = theme(
            axis_text_x=element_blank(),
            axis_ticks_major_x=element_blank(),
            strip_text_x=element_text(angle=90),
            strip_text_y=element_text(angle=0),
        )
    if occurrences:
        labels = {"labels": _key_labels}
        layout.append(scale_x_discrete(**labels) if swap_axes else scale_y_discrete(**labels))

    components = [
        geom_tile(raster=rasterized),
        *layout,
        fill_scale,
        labs(x="", y="", fill=fill_lab),
        _family_theme("matrix"),
        hide,
    ]
    threshold = annotation_threshold(annotate)
    if threshold is not None:
        style = _active_style()
        components.append(
            geom_contrast_text(
                aes(label="value"),
                min_cell_pt=threshold,
                format_string="{:.2g}",
                size=style.axis_text_size if style is not None else 7,
            )
        )
    return ggplot(tidy, aes(xcol, ycol, fill="value")) + components
