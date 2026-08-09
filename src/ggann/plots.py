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
    scale_fill_cmap,
    scale_fill_gradient2,
    scale_size,
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
from ._grouping import _downsample_cells, _group_categories, _order_groups
from ._matplotlib_backend import categorical_palette, promote_matplotlib_plot
from ._palette import obs_colors, scale_color_obs
from ._resolve import (
    ObsmRef,
    Ref,
    embedding_coords,
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
    return pd.api.types.is_numeric_dtype(series) and not isinstance(
        series.dtype, pd.CategoricalDtype
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
    backend: str = "plotnine",
    rasterized: bool = False,
    show_axes: bool | None = None,
    equal_aspect: bool | None = None,
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

    Parameters
    ----------
    adata
        Annotated data matrix.
    basis : str, default="umap"
        Embedding basis in ``adata.obsm``.
    color : str or selector, optional
        Observation column, gene, or explicit :func:`ggann.obs`/:func:`ggann.gene`.
    split_by : str, optional
        Observation column used for facets.
    layer, use_raw : optional
        Expression source for gene colour; these options are mutually exclusive.
    size, alpha : float
        Point size and opacity.
    pointdensity : bool, optional
        Use density colouring when ``color`` is absent.
    label : bool, default=False
        Label categorical colour centroids.
    label_size : float
        Centroid-label size.
    low, high : str
        End colours for continuous values.
    downsample : int, optional
        Maximum total observations to draw.
    random_state : int, optional
        Reproducible downsampling seed; ``None`` is non-deterministic.
    backend : {"plotnine", "matplotlib"}, default="plotnine"
        Rendering path. The explicit Matplotlib path accelerates an unsplit
        point scatter while returning a composable ggplot subclass.
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
        Composable embedding plot. The Matplotlib backend returns a
        :class:`ggann.MatplotlibGGPlot` subclass.

    Raises
    ------
    KeyError
        If the embedding or an explicit colour/split source is missing.
    ValueError
        If the embedding has fewer than two coordinates, downsampling/backend is
        invalid, or the Matplotlib backend is asked for an unsupported density,
        facet, or centroid-label layout.

    Notes
    -----
    Only two embedding coordinates and the requested colour metadata are projected.
    Downsampling occurs before extraction and never mutates ``adata``.

    Examples
    --------
    >>> p = plot_embedding(adata, "umap", color="cell_type")
    """
    _validate_backend(backend)
    if backend == "matplotlib" and (split_by is not None or label):
        raise ValueError(
            "backend='matplotlib' currently supports unsplit embeddings without "
            "centroid labels; use backend='plotnine' for facets or label=True."
        )
    adata = _downsample_cells(adata, None, downsample, random_state=random_state)
    key = embedding_key(adata, basis)
    width = adata.obsm[key].shape[1]
    if width < 2:
        raise ValueError(
            f"Embedding '{basis}' has only {width} dimension(s); "
            "plot_embedding requires at least 2."
        )
    if backend == "matplotlib":
        df, xcol, ycol, native_cname = _native_embedding_frame(
            adata,
            key,
            color,
            layer=layer,
            use_raw=use_raw,
        )
    else:
        coordinate_refs = [obsm(key, 0), obsm(key, 1)]
        requested = [*coordinate_refs]
        if color is not None:
            requested.append(color)
        if split_by is not None:
            requested.append(split_by)
        df = resolve_frame(adata, requested, layer=layer, use_raw=use_raw)
        xcol, ycol = (plain_name(adata, ref) for ref in coordinate_refs)
        native_cname = None
    facet = pe.facet_wrap("~" + split_by) if split_by is not None else None

    def _with_facet(components):
        return [*components, facet] if facet is not None else components

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

    if label and color is None:
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
            return ggplot(df, aes(xcol, ycol)) + _with_facet(
                [
                    pe.geom_pointdensity(size=size, alpha=alpha, raster=rasterized),
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
        if label:
            warnings.warn(
                f"plot_embedding: label=True is ignored for the numeric color {color!r} "
                "(centroid labels need a categorical color).",
                stacklevel=2,
            )
        # Draw low-expression cells first so high-expression cells are not occluded
        # (mirrors scanpy's sc.pl.embedding ordering).
        df = df.sort_values(cname)
        style = _active_style()
        finite = df[cname].to_numpy(dtype=float)
        finite = finite[np.isfinite(finite)]
        signed = bool(len(finite) and finite.min() < 0 < finite.max())
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
    components = [
        _solid_point(),
        scale_color_obs(adata, cname),
        # enlarge the legend swatches so categories stay readable (scplotter does this)
        guides(
            color=guide_legend(
                override_aes={"size": 4},
                ncol=(2 if _active_style() is not None and df[cname].nunique() > 8 else None),
                byrow=_active_style() is not None,
            )
        ),
        *_embedding_style(show_axes, equal_aspect),
    ]
    if label:
        cents = _centroid_labels(df, cname, xcol, ycol, split_by=split_by)
        # white-backed repelled labels at centroids, like scplotter's label_bg="white"
        components.append(
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
        ggplot(df, aes(xcol, ycol, color=cname)) + _with_facet(components),
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
    ncol: int | None = None,
    size: float = 1.2,
    alpha: float = 0.9,
    cmap: str = "magma",
    downsample: int | None = None,
    random_state: int | None = 0,
    rasterized: bool = False,
    show_axes: bool | None = None,
    equal_aspect: bool | None = None,
):
    """Multi-gene embedding grid: one faceted panel per feature.

    Panels share a single expression colour scale, which is best for comparing
    magnitudes *across* genes. Because the scale is shared, a low-range gene shown
    next to a high-range one will look faint -- for independent per-gene colour
    bars (like ``sc.pl.umap(color=[...])``), compose separate ``plot_embedding``
    calls with the re-exported ``Wrap`` / ``plot_layout`` instead. ``downsample=N``
    caps cells before drawing for lighter panels on large data; ``random_state``
    controls the retained cells.

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
    ncol : int, optional
        Facet columns.
    size, alpha : float
        Point size and opacity.
    cmap : str
        Matplotlib colormap name.
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
        Composable faceted feature plot.

    Raises
    ------
    KeyError
        If the embedding or an explicit source is missing.
    ValueError
        If no numeric feature resolves, the embedding is one-dimensional, or
        downsampling/source selection is invalid.

    Notes
    -----
    annplyr projects requested fields only; ``adata`` is never mutated.

    Examples
    --------
    >>> p = plot_features(adata, ["CD3D", "NKG7"], basis="umap")
    """
    features = ordered_unique(features)
    adata = _downsample_cells(adata, None, downsample, random_state=random_state)
    coords = embedding_coords(adata, basis)
    if coords.shape[1] < 2:
        raise ValueError(f"Embedding '{basis}' has fewer than 2 dimensions.")
    xcol, ycol = coords.columns[:2]

    values = resolve_frame(adata, list(features), layer=layer, use_raw=use_raw)
    feats = list(
        dict.fromkeys(f for f in features if f in values.columns and _is_numeric(values[f]))
    )
    if not feats:
        raise ValueError("plot_features needs at least one numeric feature (gene or metric).")

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
        + geom_point(size=size, alpha=alpha, raster=rasterized)
        + pe.facet_wrap("~feature", ncol=ncol)
        + _continuous_scale("color", long["expression"], cmap, signed=False)
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
    genes : sequence of str
        Genes to summarize, in display order.
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
    cmap : str
        Matplotlib colormap name.
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

    Examples
    --------
    >>> p = plot_dotplot(adata, ["CD3D", "NKG7"], group_by="cell_type")
    """
    _validate_backend(backend)
    if backend == "matplotlib" and split_by is not None:
        raise ValueError(
            "backend='matplotlib' currently supports unsplit dotplots; "
            "use backend='plotnine' for split_by facets."
        )
    genes = ordered_unique(genes)
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
    color_lab = "scaled\nexpression" if standard_scale else "mean\nexpression"
    size_kwargs = {
        "range": size_range,
        "labels": lambda xs: [f"{x:.0%}" for x in xs],
    }
    if _active_style() is not None:
        size_kwargs.update(limits=(0, 1), breaks=(0, 0.25, 0.5, 0.75, 1.0))
    components = [
        geom_point(
            aes(size="fraction", color="mean_expression"),
            raster=rasterized,
        ),
        _continuous_scale(
            "color",
            agg["mean_expression"],
            cmap,
            signed=standard_scale == "zscore",
        ),
        scale_size(**size_kwargs),
        labs(x="", y="", color=color_lab, size="fraction\nexpressing"),
        _family_theme("matrix"),
        pe.rotate_x_text(45),
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
    if split_by:
        components.append(pe.facet_wrap(f"~{split_by}"))
    plot = ggplot(agg, aes("feature", group_by)) + components
    if backend == "plotnine":
        return plot
    return promote_matplotlib_plot(
        plot,
        kind="dotplot",
        data=agg,
        x="feature",
        y=group_by,
        value="mean_expression",
        fraction="fraction",
        cmap=cmap,
        size_range=size_range,
        value_label=color_lab,
        publication_style=_active_style(),
        rasterized=rasterized,
        signed=standard_scale == "zscore",
        annotation_min_cell_pt=threshold,
    )


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
    genes : sequence of str
        Genes to summarize, in display order.
    group_by : str
        Primary observation grouping column.
    split_by : str, optional
        Additional grouping and facet column.
    layer, use_raw : optional
        Mutually exclusive expression-source selection.
    standard_scale : {None, "var", "group", "zscore"}
        Optional scaling of group means.
    cmap : str
        Matplotlib colormap name.
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
    not mutated.

    Examples
    --------
    >>> p = plot_matrixplot(adata, ["CD3D", "NKG7"], group_by="cell_type")
    """
    _validate_backend(backend)
    if backend == "matplotlib" and split_by is not None:
        raise ValueError(
            "backend='matplotlib' currently supports unsplit matrixplots; "
            "use backend='plotnine' for split_by facets."
        )
    genes = ordered_unique(genes)
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
    color_lab = "scaled\nexpression" if standard_scale else "mean\nexpression"
    components = [
        geom_tile(raster=rasterized),
        _continuous_scale(
            "fill",
            agg["mean_expression"],
            cmap,
            signed=standard_scale == "zscore",
        ),
        labs(x="", y="", fill=color_lab),
        _family_theme("matrix"),
        pe.rotate_x_text(45),
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
    if split_by:
        components.append(pe.facet_wrap(f"~{split_by}"))
    plot = ggplot(agg, aes("feature", group_by, fill="mean_expression")) + components
    if backend == "plotnine":
        return plot
    return promote_matplotlib_plot(
        plot,
        kind="matrixplot",
        data=agg,
        x="feature",
        y=group_by,
        value="mean_expression",
        cmap=cmap,
        value_label=color_lab,
        publication_style=_active_style(),
        rasterized=rasterized,
        signed=standard_scale == "zscore",
        annotation_min_cell_pt=threshold,
    )


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
        Observation column defining density facets.
    size : float
        Point size.
    ncol : int, optional
        Facet columns.
    cmap : str
        Matplotlib colormap name.
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
    by default. The input is not mutated.

    Examples
    --------
    >>> p = plot_embedding_density(adata, "umap", group_by="cell_type")
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
        plot = ggplot(df, aes(xcol, ycol, color="density")) + geom_point(
            size=size, raster=rasterized
        )
    else:
        gcol = resolve_frame(adata, [group_by])[[group_by]]
        df = coords.join(gcol)
        # groupby(...).apply concatenates in group-sorted order, so reindex back
        # to df's cell order (a bare .to_numpy() would misalign interleaved groups)
        per_group = df.groupby(group_by, observed=True, group_keys=False).apply(_density)
        df["density"] = per_group.reindex(df.index).to_numpy()
        cats = _group_categories(adata, group_by)
        df = _order_groups(df, group_by, cats)
        plot = (
            ggplot(df, aes(xcol, ycol, color="density"))
            + geom_point(size=size, raster=rasterized)
            + pe.facet_wrap("~" + group_by, ncol=ncol)
        )
    return (
        plot
        + _continuous_scale("color", df["density"], cmap, signed=False)
        + _embedding_style(show_axes, equal_aspect)
    )


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
    genes : sequence of str
        Genes shown as rows.
    group_by : str
        Observation column defining cell blocks.
    layer, use_raw : optional
        Mutually exclusive expression-source selection.
    cmap : str
        Matplotlib colormap name.
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
    adata = _downsample_cells(adata, group_by, downsample, random_state=random_state)
    genes = ordered_unique(genes)
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

    components = [
        geom_tile(raster=rasterized),
        facet_grid(". ~ " + group_by, scales="free_x", space="free_x"),
        _continuous_scale(
            "fill",
            tidy["value"],
            cmap,
            signed=standard_scale == "zscore",
        ),
        labs(x="", y="", fill=fill_lab),
        _family_theme("matrix"),
        theme(
            axis_text_x=element_blank(),
            axis_ticks_major_x=element_blank(),
            strip_text_x=element_text(angle=90),
        ),
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
    return ggplot(tidy, aes("cell_rank", "feature", fill="value")) + components
