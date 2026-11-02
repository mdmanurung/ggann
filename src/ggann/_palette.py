"""Reuse scanpy's stored category colours so ggann matches the rest of an analysis.

scanpy stores a per-category colour list in ``adata.uns[f"{col}_colors"]``, aligned
to ``adata.obs[col].cat.categories``. When present, ggann picks it up so a cluster
plotted here has the same colours as ``sc.pl.umap`` of the same object.
"""

from __future__ import annotations

import warnings
from typing import cast

import pandas as pd
import plotnine_extra as pe
from matplotlib.colors import is_color_like, to_hex
from plotnine import scale_color_manual, scale_fill_manual

from ._compat import SKIP_GGANN_FRAMES

__all__ = ["obs_colors", "scale_color_obs", "scale_colour_obs", "scale_fill_obs"]

_WARNED: set[tuple] = set()


def _warn_once(key: tuple, message: str) -> None:
    if key in _WARNED:
        return
    _WARNED.add(key)
    warnings.warn(message, UserWarning, skip_file_prefixes=SKIP_GGANN_FRAMES)


def _audit_palette(colors, source: str, *, accessibility: bool) -> None:
    """Warn about duplicate colours and, in a publication context, weak CVD separation.

    Stored palettes are never changed: the warning names the problem and the
    mapping is kept, so figures stay consistent with the rest of an analysis.
    """
    hexes = tuple(to_hex(c, keep_alpha=False) for c in colors)
    if len(set(hexes)) != len(hexes):
        _warn_once(
            ("dup", source, hexes),
            f"{source} contains duplicate colours; categories sharing a colour cannot be "
            "told apart.",
        )
    if accessibility:
        from ._palette_qa import CVD_MIN_SEPARATION, cvd_shortfalls

        short = cvd_shortfalls(hexes)
        if short:
            worst = ", ".join(f"{name} {value:.1f}" for name, value in short.items())
            _warn_once(
                ("cvd", hexes),
                f"{source}: {len(hexes)} colours fall below the colour-vision separation "
                f"threshold ({CVD_MIN_SEPARATION:g} CIEDE2000; minimum {worst}). Colours "
                "are kept as given; the first eight publication colours are the only "
                "measured guarantee, so add a redundant encoding (labels, shapes).",
            )


def _is_continuous(series: pd.Series) -> bool:
    """Numeric metadata is continuous; bool and categorical are discrete."""
    return (
        pd.api.types.is_numeric_dtype(series)
        and not pd.api.types.is_bool_dtype(series)
        and not isinstance(series.dtype, pd.CategoricalDtype)
    )


def _scanpy_default_colors(n: int) -> list[str] | None:
    """Scanpy's default categorical palette for ``n`` categories, or ``None``.

    Mirrors ``scanpy.plotting._utils.add_colors_for_categorical_sample_annotation``:
    20, 28, then 102 colours; beyond that scanpy falls back to grey, so ggann keeps
    its hue scale instead.
    """
    from scanpy.plotting import palettes

    for palette in (palettes.default_20, palettes.default_28, palettes.default_102):
        if n <= len(palette):
            return [str(c) for c in palette[:n]]
    return None


def obs_colors(adata, col: str) -> dict | None:
    """Return ``{category: hex}`` from ``adata.uns[f"{col}_colors"]``, or ``None``.

    Falls back to ``None`` when the column is not a categorical or has no stored
    colours (or the colour list is shorter than the categories, or contains an
    invalid colour). Duplicate stored colours trigger a warning and are kept.

    Parameters
    ----------
    adata
        Annotated data matrix.
    col : str
        Categorical observation column.

    Returns
    -------
    dict or None
        Category-to-colour mapping when a complete stored palette exists.

    Raises
    ------
    None
        Missing or non-categorical columns return ``None``.

    Notes
    -----
    The stored palette is read without modifying ``adata.uns``.

    Examples
    --------
    >>> palette = obs_colors(adata, "cell_type")
    """
    dtype = adata.obs[col].dtype if col in adata.obs else None
    if not isinstance(dtype, pd.CategoricalDtype):
        return None
    colors = adata.uns.get(f"{col}_colors")
    cats = list(adata.obs[col].cat.categories)
    from .publication import _active_style, publication_palette

    style = _active_style()
    source = f"adata.uns['{col}_colors']"
    if style is None:
        if colors is None or len(colors) < len(cats):
            return None
        if not all(is_color_like(str(color)) for color in colors[: len(cats)]):
            warnings.warn(
                f"{source} contains invalid colours; using the default palette.",
                UserWarning,
                skip_file_prefixes=SKIP_GGANN_FRAMES,
            )
            return None
        _audit_palette(colors[: len(cats)], source, accessibility=False)
        return {cat: str(c) for cat, c in zip(cats, colors)}

    if colors is None:
        generated = publication_palette("qualitative", categories=cats, style=style)
        _audit_palette(
            list(generated.values()), f"publication palette for '{col}'", accessibility=True
        )  # type: ignore[union-attr]
        return cast(dict[object, str], generated)
    valid = len(colors) == len(cats) and all(is_color_like(str(color)) for color in colors)
    if not valid:
        warnings.warn(
            f"adata.uns['{col}_colors'] must contain exactly one valid colour for "
            f"each of the {len(cats)} categories; using the deterministic publication palette.",
            UserWarning,
            skip_file_prefixes=SKIP_GGANN_FRAMES,
        )
        generated = publication_palette("qualitative", categories=cats, style=style)
        _audit_palette(
            list(generated.values()), f"publication palette for '{col}'", accessibility=True
        )  # type: ignore[union-attr]
        return cast(dict[object, str], generated)
    _audit_palette(colors, source, accessibility=True)
    return {cat: str(c) for cat, c in zip(cats, colors, strict=True)}


def _publication_obs_mapping(adata, col: str) -> dict | None:
    from .publication import _active_style, publication_palette

    style = _active_style()
    if style is None or col not in adata.obs:
        return None
    series = adata.obs[col]
    if isinstance(series.dtype, pd.CategoricalDtype):
        return obs_colors(adata, col)
    categories = list(pd.unique(series.dropna()))
    generated = publication_palette("qualitative", categories=categories, style=style)
    _audit_palette(
        list(generated.values()),  # type: ignore[union-attr]
        f"publication palette for '{col}'",
        accessibility=True,
    )
    return cast(dict[object, str], generated)


def _obs_mapping(adata, series: pd.Series, col: str) -> dict | None:
    """Stored palette, else the publication palette, else scanpy's default palette."""
    mapping = obs_colors(adata, col) or _publication_obs_mapping(adata, col)
    if mapping is not None:
        return mapping
    categories = list(pd.Categorical(series.dropna()).categories)
    colors = _scanpy_default_colors(len(categories))
    return None if colors is None else dict(zip(categories, colors))


def annotation_colors(adata, ann_df: pd.DataFrame) -> dict[str, dict] | None:
    """Per-column category colours for ``HeatmapAnnotation(colors=...)``, or ``None``.

    PyComplexHeatmap needs a mapping for *every* annotation column or none, so
    this returns ``None`` unless every column is a discrete ``obs`` column with a
    stored, publication or scanpy-default palette. ``adata.uns`` is not modified.
    """
    result = {}
    for col in ann_df.columns:
        if col not in adata.obs or _is_continuous(adata.obs[col]):
            return None
        mapping = _obs_mapping(adata, adata.obs[col], col)
        if mapping is None:
            return None
        values = ann_df[col].dropna().unique()
        if any(value not in mapping for value in values):
            return None
        result[col] = {value: mapping[value] for value in values}
    return result or None


def _obs_scale(aesthetic: str, adata, col: str, kwargs: dict):
    """Shared dispatch for the colour/fill observation scales (validated here)."""
    if col not in adata.obs:
        raise KeyError(f"obs column {col!r} not found; cannot build a {aesthetic} scale.")
    series = adata.obs[col]
    from .publication import _active_style

    style = _active_style()
    if _is_continuous(series):
        if style:
            kwargs.setdefault("na_value", style.missing_color)
        kwargs.setdefault("cmap_name", style.sequential_cmap if style else "viridis")
        return (pe.scale_color_cmap if aesthetic == "color" else pe.scale_fill_cmap)(**kwargs)
    mapping = _obs_mapping(adata, series, col)
    if mapping is None:
        return (pe.scale_color_hue if aesthetic == "color" else pe.scale_fill_hue)(**kwargs)
    if style:
        kwargs.setdefault("na_value", style.missing_color)
        kwargs.setdefault("drop", False)
    manual = scale_color_manual if aesthetic == "color" else scale_fill_manual
    return manual(values=mapping, **kwargs)


def scale_color_obs(adata, col: str, **kwargs):
    """Build a categorical colour scale from observation metadata.

    Parameters
    ----------
    adata
        Annotated data matrix.
    col : str
        Observation column.
    **kwargs
        Passed to the plotnine manual, hue or colour-map scale.

    Returns
    -------
    plotnine.scales.scale
        Composable colour scale.

    Raises
    ------
    KeyError
        If ``col`` is not an observation column.
    ValueError
        If plotnine rejects a forwarded scale argument.

    Notes
    -----
    Numeric columns get a continuous colour-map scale (viridis, or the active
    publication style's sequential map); boolean and categorical columns get a
    discrete scale. Discrete colours come from the stored
    ``adata.uns["<col>_colors"]`` when complete, then the publication palette, then
    scanpy's default palette (up to 102 categories, hues beyond). The palette is
    assigned in category order, so reordering categories reassigns colours.
    ``adata.uns`` is never modified, and a scale added later replaces this one.

    Examples
    --------
    >>> p = p + scale_color_obs(adata, "cell_type")
    """
    return _obs_scale("color", adata, col, kwargs)


# British-spelling alias, mirroring plotnine.
scale_colour_obs = scale_color_obs


def scale_fill_obs(adata, col: str, **kwargs):
    """Build a categorical fill scale from observation metadata.

    Parameters
    ----------
    adata
        Annotated data matrix.
    col : str
        Observation column.
    **kwargs
        Passed to the plotnine manual, hue or colour-map scale.

    Returns
    -------
    plotnine.scales.scale
        Composable fill scale.

    Raises
    ------
    KeyError
        If ``col`` is not an observation column.
    ValueError
        If plotnine rejects a forwarded scale argument.

    Notes
    -----
    Numeric columns get a continuous colour-map scale (viridis, or the active
    publication style's sequential map); boolean and categorical columns get a
    discrete scale. Discrete colours come from the stored
    ``adata.uns["<col>_colors"]`` when complete, then the publication palette, then
    scanpy's default palette (up to 102 categories, hues beyond). The palette is
    assigned in category order, so reordering categories reassigns colours.
    ``adata.uns`` is never modified, and a scale added later replaces this one.

    Examples
    --------
    >>> p = p + scale_fill_obs(adata, "cell_type")
    """
    return _obs_scale("fill", adata, col, kwargs)
