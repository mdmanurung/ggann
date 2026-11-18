"""Radar layers: an ordered, closed polygon per series on a radial coord."""

from __future__ import annotations

import numpy as np
import pandas as pd
from mizani.breaks import breaks_extended
from plotnine import (
    aes,
    coord_radial,
    element_blank,
    element_line,
    geom_point,
    geom_polygon,
    ggplot,
    labs,
    scale_x_discrete,
    scale_y_continuous,
    theme,
)

from .._plotnine_compat import stat
from ..publication import _family_theme
from ._common import palette_scales, prepare_levels

__all__ = ["coord_radar", "geom_radar", "ggradar", "stat_radar"]


class stat_radar(stat):
    """Sort each series by axis and close it, so a polygon or path is whole.

    Parameters
    ----------
    mapping : plotnine.mapping.aes, optional
        Aesthetics. ``x`` (the axis, discrete) and ``y`` (the radius) are
        required. Each ``group`` is one series.
    data : pandas.DataFrame, optional
        Layer data. Defaults to the plot data.
    **kwargs
        ``closed`` (bool, default True) appends the first row of each series
        at its end. ``geom`` is the geom that draws the series (default
        ``geom_polygon``). ``na_rm`` and ``position`` work as in any plotnine
        stat.

    Returns
    -------
    stat_radar
        A layer to add to a ``plotnine.ggplot``. Its output is the input sorted
        by ``x`` within each group, plus the first row when ``closed=True``.

    Raises
    ------
    plotnine.exceptions.PlotnineError
        If a keyword matches no stat, geom or aesthetic parameter.

    Notes
    -----
    The layer tier draws what it gets. A series with a missing axis has one
    vertex fewer and no warning is raised. Draw it with ``coord_radar``, a
    discrete x scale with ``scale_x_discrete(expand=(0, 0.5))`` so the axes are
    evenly spaced, and r limits that start at 0. For K axes the recipe is
    ``geom_radar() + scale_x_discrete(expand=(0, 0.5)) +
    coord_radar(start=-pi / K)``; ``start=-pi / K`` puts the first axis at 12
    o'clock. Values outside the r limits are censored and their vertices
    vanish. Once ggann is imported, ``stat="radar"`` in any plotnine layer
    resolves to this class.

    Examples
    --------
    >>> import pandas as pd
    >>> from plotnine import aes, ggplot
    >>> import ggann as ag
    >>> df = pd.DataFrame({"axis": list("abcde"), "v": [4, 2, 3, 1, 2], "s": "s1"})
    >>> p = ggplot(df, aes("axis", "v", group="s")) + ag.stat_radar()
    """

    REQUIRED_AES = {"x", "y"}
    DEFAULT_PARAMS = {"geom": geom_polygon, "closed": True}

    def compute_group(self, data, scales):
        out = data.sort_values("x", kind="stable")
        if self.params["closed"] and len(out):
            out = pd.concat([out, out.iloc[[0]]], ignore_index=True)
        return out


class geom_radar(geom_polygon):
    """Draw each series as a translucent polygon, with ``stat_radar`` as the stat.

    Parameters
    ----------
    mapping : plotnine.mapping.aes, optional
        Aesthetics. ``x`` (the axis) and ``y`` (the radius) are required;
        ``group``, ``fill`` and ``color`` mark the series.
    data : pandas.DataFrame, optional
        Layer data. Defaults to the plot data.
    **kwargs
        ``closed`` goes to ``stat_radar``. Other keywords are polygon
        aesthetics or parameters. ``alpha`` defaults to 0.25.

    Returns
    -------
    geom_radar
        A layer to add to a ``plotnine.ggplot``.

    Raises
    ------
    plotnine.exceptions.PlotnineError
        If a keyword matches no stat, geom or aesthetic parameter.

    Notes
    -----
    Edges stay straight only under ``coord_radar``; ``coord_radial`` draws each
    edge as a spiral. See ``stat_radar`` for the recipe and the limits.

    Examples
    --------
    >>> import pandas as pd
    >>> from plotnine import aes, ggplot, scale_x_discrete
    >>> import ggann as ag
    >>> df = pd.DataFrame({"axis": list("abcde"), "v": [4, 2, 3, 1, 2], "s": "s1"})
    >>> p = (
    ...     ggplot(df, aes("axis", "v", group="s", fill="s"))
    ...     + ag.geom_radar()
    ...     + scale_x_discrete(expand=(0, 0.5))
    ...     + ag.coord_radar(start=-0.2 * 3.14159)
    ... )
    """

    DEFAULT_PARAMS = {"stat": stat_radar}
    DEFAULT_AES = {**geom_polygon.DEFAULT_AES, "alpha": 0.25}


class coord_radar(coord_radial):
    """A radial coord that keeps straight edges straight.

    ``coord_radial`` draws each polygon edge as a spiral. This subclass turns
    that off for every geom in the plot, so a radar polygon has one vertex per
    axis.

    Parameters
    ----------
    theta : str
        ``"x"`` or ``"y"``: the aesthetic mapped to the angle.
    start : float
        Angle of the start of the theta range, in radians, clockwise from 12
        o'clock.
    end : float, optional
        Angle of the end of the theta range. ``None`` is a full circle.
    expand : bool
        Whether the scales add their usual padding.
    inner_radius : float
        Hole radius as a fraction of the outer radius.
    rotate_angle : bool
        Whether text follows the angle.
    thetalim : tuple of float, optional
        Limits of the theta scale.
    rlim : tuple of float, optional
        Limits of the radius scale.
    reverse : {"none", "theta", "r", "thetar"}
        Which direction to reverse.

    Returns
    -------
    coord_radar
        A coord to add to a ``plotnine.ggplot``.

    Raises
    ------
    ValueError
        If ``reverse`` is not one of the options.
    plotnine.exceptions.PlotnineError
        If ``start`` and ``end`` are equal.

    Notes
    -----
    Use it for radar plots only. Pie wedges and any curved path under this coord
    draw as straight chords. It relies on ``coord_radial.is_linear`` and
    ``coord_radial.munch``, which a canary test watches.

    Examples
    --------
    >>> import ggann as ag
    >>> coord = ag.coord_radar(start=-0.6283)
    """

    is_linear = True

    def munch(self, data, panel_params):
        return data


def _radar_plot(
    frame: pd.DataFrame,
    axis: str,
    value: str,
    series: str | None,
    n_axes: int,
    rlim: tuple[float, float],
    fill_alpha: float,
    points: bool,
):
    """Layers, scales, coord and theme shared by ``ggradar`` and ``plot_proportions``."""
    if series is None:
        mapping = aes(axis, value, group=1)
        style = {"color": "#444444", "fill": "#444444"}
    else:
        mapping = aes(axis, value, group=series, color=series, fill=series)
        style = {}
    plot = ggplot(frame, mapping) + geom_radar(alpha=fill_alpha, **style)
    if points:
        plot += geom_point(size=1.5, **{"color": "#444444"} if series is None else {})
    # r must start at 0, and a break at the lower limit shifts the tick labels
    # by one, so keep only the breaks above it
    return (
        plot
        + scale_x_discrete(expand=(0, 0.5))
        + scale_y_continuous(
            limits=rlim,
            expand=(0, 0),
            breaks=lambda lim: [b for b in breaks_extended()(lim) if b > lim[0]],
        )
        + coord_radar(start=-np.pi / n_axes)
        + labs(color=series, fill=series)
        + _family_theme("standard")
        + theme(
            panel_border=element_blank(),
            axis_title=element_blank(),
            panel_grid_major=element_line(color="#d9d9d9", size=0.3),
        )
    )


def ggradar(
    df: pd.DataFrame,
    axis: str,
    value: str,
    series: str | None = None,
    *,
    axis_order=None,
    series_order=None,
    rlim: tuple[float, float] | None = None,
    fill_alpha: float = 0.25,
    points: bool = True,
    palette=None,
):
    """Draw a radar plot: one polygon per series, one spoke per axis level.

    Parameters
    ----------
    df : pandas.DataFrame
        One row per ``(axis, series)`` pair.
    axis : str
        Column whose levels are the spokes. At least 3 levels are needed. The
        first level sits at 12 o'clock and the rest follow clockwise.
    value : str
        Numeric column drawn as the radius.
    series : str, optional
        Column whose levels are the polygons. ``None`` draws one polygon.
    axis_order, series_order : sequence, optional
        Level orders. ``None`` uses the categorical order, then sorted values.
    rlim : tuple of float, optional
        Radius limits ``(low, high)``. ``None`` gives ``(0, max(value))``.
    fill_alpha : float
        Opacity of the polygon fill.
    points : bool
        Whether to draw a point at each vertex.
    palette : dict, list or str, optional
        Colours for the series. A dict must cover every level. ``None`` keeps
        plotnine's default.

    Returns
    -------
    plotnine.ggplot
        A plot that accepts further layers, scales and themes.

    Raises
    ------
    KeyError
        If a named column is not in ``df``.
    ValueError
        If there are fewer than 3 axis levels, a pair is duplicated or missing,
        ``axis`` or ``series`` has a missing label, a value is missing or not
        finite, ``rlim`` is not an increasing pair, a value lies outside
        ``rlim``, or the palette is incomplete.

    Notes
    -----
    The radius is the value itself: there is no rescaling per axis, so a
    spoke with small values stays near the centre. The default limits start at
    0, and so reject negative values. A value outside ``rlim`` raises, because
    the plot would censor it and the polygon would lose a vertex. Spokes are
    evenly spaced and the grid is circular.

    Examples
    --------
    >>> import pandas as pd
    >>> import ggann as ag
    >>> df = pd.DataFrame({"axis": list("abcabc"), "v": [1, 2, 3, 2, 1, 2], "s": list("xxxyyy")})
    >>> p = ag.ggradar(df, "axis", "v", series="s")
    """
    keep = [c for c in (axis, value, series) if c is not None]
    frame, axis_levels = prepare_levels(df, axis, axis_order, keep=keep)
    if len(axis_levels) < 3:
        raise ValueError(f"a radar plot needs at least 3 axis levels, got {len(axis_levels)}")
    if series is not None:
        frame, _ = prepare_levels(frame, series, series_order, keep=keep)
    labels = [c for c in (axis, series) if c is not None]
    if bool(frame[labels].isna().to_numpy().any()):
        raise ValueError(f"{labels} must have no missing labels")
    frame[value] = pd.to_numeric(frame[value], errors="raise").astype(float)
    if not np.isfinite(frame[value]).all():
        raise ValueError(f"{value!r} has a missing or non-finite value; remove those rows first")
    pairs = labels
    if frame.duplicated(pairs).any():
        raise ValueError(f"more than one row for the same {pairs}; give one row per pair")
    series_levels = [None] if series is None else list(frame[series].cat.categories)
    if len(frame) != len(axis_levels) * len(series_levels):
        raise ValueError(f"every {pairs} pair needs a row; some are missing")
    if rlim is None:
        rlim = (0.0, float(frame[value].max()))
    rlim = (float(rlim[0]), float(rlim[1]))
    if not rlim[0] < rlim[1]:
        raise ValueError(f"rlim must be an increasing pair, got {rlim!r}")
    outside = (frame[value] < rlim[0]) | (frame[value] > rlim[1])
    if bool(outside.any()):
        raise ValueError(
            f"{int(outside.sum())} value(s) of {value!r} lie outside rlim {rlim}; "
            "pass a wider rlim or remove them"
        )
    plot = _radar_plot(frame, axis, value, series, len(axis_levels), rlim, fill_alpha, points)
    if series is None:
        return plot
    return plot + palette_scales(palette, series_levels, ("color", "fill"))
