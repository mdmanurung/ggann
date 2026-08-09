"""Default theme and single-cell colour scales for ggann.

These are thin conveniences over plotnine / plotnine-extra so that ggann figures
have a consistent, publication-friendly look while remaining ordinary plotnine
objects the user can override with any ``+ theme(...)`` / ``+ scale_*``.
"""

from __future__ import annotations

import plotnine_extra as pe
from plotnine import (
    element_blank,
    element_line,
    element_rect,
    element_text,
    theme,
    theme_bw,
    theme_gray,
    theme_set,
)

__all__ = [
    "theme_ggann",
    "set_theme",
    "reset_theme",
    "sizes",
    "scale_color_expression",
    "scale_colour_expression",
    "scale_fill_expression",
    "scale_color_celltype",
    "scale_colour_celltype",
    "scale_fill_celltype",
]

# ggplot2 points-to-mm constant: geom_text `size` is in mm, theme text in pt.
_PT_TO_MM = 1 / 2.845276


class _Sizes:
    """A font-size scale (in **pt**) kept in sync with :func:`set_theme`.

    Mirrors exactplot's ``xp$fontsize``: use these so text in ``geom_text`` /
    annotations matches the theme. Because plotnine's ``geom_text`` ``size`` is in
    mm (a ggplot2 footgun), convert with :meth:`geom`::

        geom_text(size=ag.sizes.geom(ag.sizes.small))

    Attributes
    ----------
    normal, small, large, title : float
        Related text sizes in points.

    Notes
    -----
    :func:`set_theme` updates this process-global scale. ``geom`` converts points
    to the millimetre unit used by plotnine text geometries.

    Examples
    --------
    >>> sizes.geom(sizes.small) > 0
    True
    """

    def __init__(self, base_size: float = 11) -> None:
        self.update(base_size)

    def update(self, base_size: float) -> None:
        self.normal = float(base_size)
        self.small = round(base_size * 0.8, 2)
        self.large = round(base_size * 1.25, 2)
        self.title = round(base_size * 1.3, 2)

    @staticmethod
    def geom(pt: float) -> float:
        """Convert a point size to plotnine's ``geom_text`` size unit (mm)."""
        return pt * _PT_TO_MM

    def __repr__(self) -> str:
        return (
            f"Sizes(normal={self.normal}, small={self.small}, "
            f"large={self.large}, title={self.title})"
        )


#: Font-size scale, refreshed by :func:`set_theme` (see :class:`_Sizes`).
sizes = _Sizes()


def theme_ggann(base_size: float = 11, base_family: str | None = None):
    """A clean, boxed theme emulating scplotter / plotthis ``theme_this``.

    The signature scplotter look is a thin black panel border (a full box) with
    no axis lines and no grid -- not the L-shaped axes of ``theme_classic``.

    Parameters
    ----------
    base_size : float
        Base text size in points.
    base_family : str, optional
        Installed font family.

    Returns
    -------
    plotnine.theme
        Composable theme object.

    Raises
    ------
    ValueError
        If plotnine rejects a theme argument.

    Notes
    -----
    Constructing the theme does not change global state.

    Examples
    --------
    >>> p = p + theme_ggann(base_size=9)
    """
    result = theme_bw(base_size=base_size, base_family=base_family)
    # ``result`` is fresh, so in-place composition avoids an otherwise wasted
    # deep copy of the complete base theme on every helper construction.
    result += theme(
        panel_grid=element_blank(),
        panel_border=element_rect(color="black", fill=None, size=1),
        axis_line=element_blank(),
        axis_ticks=element_line(color="black"),
        axis_text=element_text(color="black"),
        legend_key=element_blank(),
        strip_background=element_blank(),
        strip_text=element_text(color="black"),
    )
    return result


def set_theme(base_size: float = 11, family: str | None = None, *, register: bool = True):
    """Make :func:`theme_ggann` the default look for every figure (exactplot ``xp_init``).

    Sets ``ggann.sizes`` from ``base_size`` and, when ``register=True`` (default),
    registers the theme as plotnine's **global default** via ``theme_set`` -- so a
    bare ``ggplot(...)`` / ``gganndata(...)`` gets the ggann look without adding
    ``+ theme_ggann()``. The theme is also returned, so it stays composable.

    This mutates plotnine's global default theme (a deliberate, process-wide side
    effect); call :func:`reset_theme` to restore plotnine's default, or pass
    ``register=False`` to only get the theme and update ``sizes``.

    ``family`` sets the font (any installed family, e.g. ``"Arial"`` /
    ``"IBM Plex Sans"``); ``None`` keeps matplotlib's default -- ggann requires no
    particular font to be installed.

    Parameters
    ----------
    base_size : float
        Base text size in points.
    family : str, optional
        Installed font family.
    register : bool, default=True
        Register the theme as plotnine's global default.

    Returns
    -------
    plotnine.theme
        Configured composable theme.

    Raises
    ------
    ValueError
        If plotnine rejects a theme argument.

    Notes
    -----
    Always updates :data:`sizes`; with ``register=True`` it deliberately mutates
    process-global plotnine theme state, but never AnnData.

    Examples
    --------
    >>> th = set_theme(base_size=9, register=False)
    """
    th = theme_ggann(base_size=base_size, base_family=family)
    sizes.update(base_size)
    if register:
        theme_set(th)
    return th


def reset_theme() -> None:
    """Restore plotnine's default theme.

    Parameters
    ----------
    None

    Returns
    -------
    None

    Raises
    ------
    RuntimeError
        If plotnine cannot replace its global theme.

    Notes
    -----
    This intentionally changes process-global plotnine state.

    Examples
    --------
    >>> reset_theme()
    """
    theme_set(theme_gray())


def scale_color_expression(cmap: str = "Reds", **kwargs):
    """Build a continuous colour scale for expression.

    Parameters
    ----------
    cmap : str
        Matplotlib colormap name.
    **kwargs
        Passed to plotnine's colour-map scale.

    Returns
    -------
    plotnine.scales.scale
        Composable colour scale.

    Raises
    ------
    ValueError
        If the colormap or a scale argument is invalid.

    Notes
    -----
    ``scale_colour_expression`` is an exact British-spelling alias.

    Examples
    --------
    >>> p = p + scale_color_expression("magma")
    """
    # Import lazily because publication styling itself builds on ``theme_ggann``.
    # A user scale added after this one still wins through plotnine's grammar.
    from .publication import _active_style

    style = _active_style()
    if style is not None:
        if cmap == "Reds":
            cmap = style.sequential_cmap
        kwargs.setdefault("na_value", style.missing_color)
    return pe.scale_color_cmap(cmap_name=cmap, **kwargs)


# British-spelling alias, mirroring plotnine.
scale_colour_expression = scale_color_expression


def scale_fill_expression(cmap: str = "viridis", **kwargs):
    """Build a continuous fill scale for expression heatmaps.

    Parameters
    ----------
    cmap : str
        Matplotlib colormap name.
    **kwargs
        Passed to plotnine's fill-map scale.

    Returns
    -------
    plotnine.scales.scale
        Composable fill scale.

    Raises
    ------
    ValueError
        If the colormap or a scale argument is invalid.

    Notes
    -----
    The scale does not inspect or copy AnnData.

    Examples
    --------
    >>> p = p + scale_fill_expression("viridis")
    """
    from .publication import _active_style

    style = _active_style()
    if style is not None:
        if cmap == "viridis":
            cmap = style.sequential_cmap
        kwargs.setdefault("na_value", style.missing_color)
    return pe.scale_fill_cmap(cmap_name=cmap, **kwargs)


def scale_color_celltype(**kwargs):
    """Build a categorical colour scale for cell identities.

    Parameters
    ----------
    **kwargs
        Passed to plotnine's hue scale.

    Returns
    -------
    plotnine.scales.scale
        Composable colour scale.

    Raises
    ------
    ValueError
        If a scale argument is invalid.

    Notes
    -----
    ``scale_colour_celltype`` is an exact British-spelling alias.

    Examples
    --------
    >>> p = p + scale_color_celltype()
    """
    return pe.scale_color_hue(**kwargs)


# British-spelling alias, mirroring plotnine.
scale_colour_celltype = scale_color_celltype


def scale_fill_celltype(**kwargs):
    """Build a categorical fill scale for cell identities.

    Parameters
    ----------
    **kwargs
        Passed to plotnine's hue scale.

    Returns
    -------
    plotnine.scales.scale
        Composable fill scale.

    Raises
    ------
    ValueError
        If a scale argument is invalid.

    Notes
    -----
    The scale does not inspect or copy AnnData.

    Examples
    --------
    >>> p = p + scale_fill_celltype()
    """
    return pe.scale_fill_hue(**kwargs)
