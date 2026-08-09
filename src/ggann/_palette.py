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
from matplotlib.colors import is_color_like
from plotnine import scale_color_manual, scale_fill_manual

__all__ = ["obs_colors", "scale_color_obs", "scale_colour_obs", "scale_fill_obs"]


def obs_colors(adata, col: str) -> dict | None:
    """Return ``{category: hex}`` from ``adata.uns[f"{col}_colors"]``, or ``None``.

    Falls back to ``None`` when the column is not a categorical or has no stored
    colours (or the colour list is shorter than the categories).

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
    if style is None:
        if colors is None or len(colors) < len(cats):
            return None
        return {cat: str(c) for cat, c in zip(cats, colors)}

    if colors is None:
        return cast(
            dict[object, str],
            publication_palette("qualitative", categories=cats, style=style),
        )
    valid = len(colors) == len(cats) and all(is_color_like(str(color)) for color in colors)
    if not valid:
        warnings.warn(
            f"adata.uns['{col}_colors'] must contain exactly one valid colour for "
            f"each of the {len(cats)} categories; using the deterministic publication palette.",
            UserWarning,
            stacklevel=2,
        )
        return cast(
            dict[object, str],
            publication_palette("qualitative", categories=cats, style=style),
        )
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
    return cast(
        dict[object, str],
        publication_palette("qualitative", categories=categories, style=style),
    )


def scale_color_obs(adata, col: str, **kwargs):
    """Build a categorical colour scale from observation metadata.

    Parameters
    ----------
    adata
        Annotated data matrix.
    col : str
        Categorical observation column.
    **kwargs
        Passed to the plotnine manual or hue scale.

    Returns
    -------
    plotnine.scales.scale
        Composable colour scale.

    Raises
    ------
    ValueError
        If plotnine rejects a forwarded scale argument.

    Notes
    -----
    Scanpy colours are reused when complete; otherwise plotnine chooses hues.

    Examples
    --------
    >>> p = p + scale_color_obs(adata, "cell_type")
    """
    mapping = obs_colors(adata, col) or _publication_obs_mapping(adata, col)
    if mapping is None:
        return pe.scale_color_hue(**kwargs)
    from .publication import _active_style

    if style := _active_style():
        kwargs.setdefault("na_value", style.missing_color)
        kwargs.setdefault("drop", False)
    return scale_color_manual(values=mapping, **kwargs)


# British-spelling alias, mirroring plotnine.
scale_colour_obs = scale_color_obs


def scale_fill_obs(adata, col: str, **kwargs):
    """Build a categorical fill scale from observation metadata.

    Parameters
    ----------
    adata
        Annotated data matrix.
    col : str
        Categorical observation column.
    **kwargs
        Passed to the plotnine manual or hue scale.

    Returns
    -------
    plotnine.scales.scale
        Composable fill scale.

    Raises
    ------
    ValueError
        If plotnine rejects a forwarded scale argument.

    Notes
    -----
    Scanpy colours are reused when complete; otherwise plotnine chooses hues.

    Examples
    --------
    >>> p = p + scale_fill_obs(adata, "cell_type")
    """
    mapping = obs_colors(adata, col) or _publication_obs_mapping(adata, col)
    if mapping is None:
        return pe.scale_fill_hue(**kwargs)
    from .publication import _active_style

    if style := _active_style():
        kwargs.setdefault("na_value", style.missing_color)
        kwargs.setdefault("drop", False)
    return scale_fill_manual(values=mapping, **kwargs)
