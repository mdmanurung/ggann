"""Pie and ring layers: one stat that turns weights into radial bands."""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
from plotnine import (
    aes,
    after_stat,
    coord_radial,
    element_blank,
    facet_wrap,
    geom_rect,
    geom_text,
    ggplot,
    labs,
    theme,
)

from .._compat import SKIP_GGANN_FRAMES, accepts_scanpy_names
from .._plotnine_compat import stat
from ..publication import _family_theme
from ._common import drop_missing_rows, palette_scales, prepare_levels

__all__ = ["geom_pie", "ggpie", "ggring", "stat_pie"]

_COLUMNS = ["xmin", "xmax", "ymin", "ymax", "y", "x", "count", "prop"]


class stat_pie(stat):
    """Turn weights into radial bands: one row per ``(x, group)`` slice.

    Under ``coord_radial(theta="y", expand=False)`` slice ``g`` spans angles
    ``2 * pi * ymin`` to ``2 * pi * ymax``, clockwise from 12 o'clock. Each
    ``x`` value is one ring, and every ring sums to 1. A ring is a pie with
    ``inner_radius > 0`` on the coord.

    Parameters
    ----------
    mapping : plotnine.mapping.aes, optional
        Aesthetics. All are optional: ``weight`` (default 1 per row), ``x``
        (default 1, one ring), ``fill`` and ``group``.
    data : pandas.DataFrame, optional
        Layer data. Defaults to the plot data.
    **kwargs
        ``width`` (float, default 1.0) is the radial width of a ring in ``x``
        units. ``label_position`` (float, default 0.5) places the output ``x``
        inside the band, as a fraction of ``width`` from the inner edge.
        ``drop_zero`` (bool, default False) removes zero-weight slices.
        ``geom`` is the geom that draws the bands (default ``geom_rect``).
        ``na_rm`` and ``position`` work as in any plotnine stat.

    Returns
    -------
    stat_pie
        A layer to add to a ``plotnine.ggplot``. Its output has ``xmin``,
        ``xmax``, ``ymin``, ``ymax``, ``y`` (the mid fraction), ``x`` (the label
        radius), ``count`` (summed weight), ``prop`` (``count`` over the ring
        total), ``group`` and ``PANEL``. Input columns that are constant within
        a slice are kept.

    Raises
    ------
    ValueError
        When the plot is built, if a weight is negative.
    plotnine.exceptions.PlotnineError
        If a keyword matches no stat, geom or aesthetic parameter.

    Notes
    -----
    The slice key is ``group``, so a continuous ``fill`` gives one slice with
    ``group = -1``. Rows with a missing weight are removed with a warning
    (``na_rm=True`` silences it). A ring with total weight 0 is dropped with a
    ``UserWarning``. Zero-weight slices stay as zero-width rows unless
    ``drop_zero=True``. Put labels on a text layer with ``drop_zero=True``;
    otherwise a zero slice gets a ``0%`` label on a slice edge.

    The default ``coord_radial(theta="y")`` has ``expand=True``, so the wedges
    cover 0.909 of the circle. Use ``coord_radial(theta="y", expand=False)``.
    Under ``coord_radar`` a wedge is drawn with straight chords; ``coord_radar``
    is for radar plots only. Once ggann is imported, ``stat="pie"`` in any
    plotnine layer resolves to this class.

    Examples
    --------
    >>> import pandas as pd
    >>> from plotnine import aes, coord_radial, ggplot
    >>> import ggann as ag
    >>> df = pd.DataFrame({"kind": list("abcd"), "n": [4, 3, 2, 1]})
    >>> p = (
    ...     ggplot(df, aes(fill="kind", weight="n"))
    ...     + ag.stat_pie()
    ...     + coord_radial(theta="y", expand=False)
    ... )
    """

    DEFAULT_AES = {"weight": 1, "x": 1}
    NON_MISSING_AES = {"weight"}
    DEFAULT_PARAMS = {
        "geom": geom_rect,
        "position": "identity",
        "na_rm": False,
        "width": 1.0,
        "label_position": 0.5,
        "drop_zero": False,
    }
    DROPPED_AES = ["weight"]
    CREATES = set(_COLUMNS)

    def compute_panel(self, data, scales):
        width = float(self.params["width"])
        label_position = float(self.params["label_position"])
        if bool((data["weight"] < 0).any()):
            raise ValueError("weights must be non-negative; found a negative weight")
        keys = ["x", "group"]
        rest = [c for c in data.columns if c not in {*keys, "weight"}]
        grouped = data.groupby(keys, sort=True, observed=True)
        # keep only the columns that are constant inside every slice
        constant = [c for c in rest if bool((grouped[c].nunique(dropna=False) <= 1).all())]
        slices = grouped[constant].first().reset_index()
        slices["count"] = grouped["weight"].sum().to_numpy(dtype=float)

        rings = []
        for x, ring in slices.groupby("x", sort=True):
            total = float(ring["count"].sum())
            if total <= 0:
                warnings.warn(
                    f"stat_pie: the ring at x = {x} has total weight 0 and is dropped",
                    UserWarning,
                    skip_file_prefixes=SKIP_GGANN_FRAMES,
                )
                continue
            ring = ring.copy()
            # cumulate counts and divide once, so the last ymax is exactly 1
            ring["ymax"] = ring["count"].cumsum().to_numpy() / total
            ring["ymin"] = ring["ymax"].shift(fill_value=0.0)
            ring["prop"] = ring["count"] / total
            rings.append(ring)
        if not rings:
            return slices.iloc[:0].assign(**{c: pd.Series(dtype=float) for c in _COLUMNS})
        out = pd.concat(rings, ignore_index=True)
        out["y"] = (out["ymin"] + out["ymax"]) / 2
        out["xmin"] = out["x"] - width / 2
        out["xmax"] = out["x"] + width / 2
        out["x"] = out["xmin"] + label_position * width
        if self.params["drop_zero"]:
            out = out[out["count"] > 0].reset_index(drop=True)
        return out


class geom_pie(geom_rect):
    """Draw pie slices or rings from weights, with ``stat_pie`` as the stat.

    A rectangle geom whose default stat computes the slice bounds. Add
    ``coord_radial(theta="y", expand=False)`` to turn the bands into a pie, and
    ``inner_radius`` on the coord to open a hole.

    Parameters
    ----------
    mapping : plotnine.mapping.aes, optional
        Aesthetics. ``fill``, ``weight``, ``x`` and ``group`` are used; none is
        required.
    data : pandas.DataFrame, optional
        Layer data. Defaults to the plot data.
    **kwargs
        ``width``, ``label_position`` and ``drop_zero`` go to ``stat_pie``.
        Other keywords are rectangle aesthetics or parameters, such as
        ``color`` and ``size`` for the outline.

    Returns
    -------
    geom_pie
        A layer to add to a ``plotnine.ggplot``.

    Raises
    ------
    ValueError
        When the plot is built, if a weight is negative.
    plotnine.exceptions.PlotnineError
        If a keyword matches no stat, geom or aesthetic parameter.

    Notes
    -----
    The layer-tier recipe is ``geom_pie() + coord_radial(theta="y",
    expand=False)``. The default coord expands the range, so the wedges would
    cover 0.909 of the circle. Do not use ``coord_radar`` with this layer.
    Several ``x`` levels draw concentric rings. See ``stat_pie`` for the
    computed columns and for missing and zero weights.

    Examples
    --------
    >>> import pandas as pd
    >>> from plotnine import aes, coord_radial, ggplot
    >>> import ggann as ag
    >>> df = pd.DataFrame({"kind": list("abcd"), "n": [4, 3, 2, 1]})
    >>> p = (
    ...     ggplot(df, aes(fill="kind", weight="n"))
    ...     + ag.geom_pie(color="white")
    ...     + coord_radial(theta="y", expand=False, inner_radius=0.5)
    ... )
    """

    DEFAULT_PARAMS = {"stat": stat_pie}


_LABELS = ("none", "count", "percent")


@accepts_scanpy_names
def ggpie(
    df: pd.DataFrame,
    group_by: str,
    value: str | None = None,
    *,
    split_by: str | None = None,
    order=None,
    label: str = "none",
    inner_radius: float = 0.0,
    start: float = 0.0,
    palette=None,
    border_color: str = "white",
    border_size: float = 0.15,
):
    """Draw a pie of the shares of ``group_by`` in a data frame.

    Parameters
    ----------
    df : pandas.DataFrame
        One row per observation, or one row per category when ``value`` is
        given.
    group_by : str
        Column whose levels are the slices.
    value : str, optional
        Column of non-negative weights. ``None`` counts rows.
    split_by : str, optional
        Column that splits the data into one pie per level (facets).
    order : sequence, optional
        Slice order. ``None`` uses the categorical order of ``group_by``, then
        sorted values.
    label : {"none", "count", "percent"}
        Text drawn in each slice. Zero slices get no label.
    inner_radius : float
        Hole radius as a fraction of the outer radius, in ``[0, 1)``.
    start : float
        Angle of the first slice edge in radians, clockwise from 12 o'clock.
    palette : dict, list or str, optional
        A dict must cover every level. A list or a Matplotlib colormap name is
        used in level order. ``None`` keeps plotnine's default.
    border_color : str
        Colour of the slice outlines.
    border_size : float
        Width of the slice outlines.

    Returns
    -------
    plotnine.ggplot
        A plot that accepts further layers, scales and themes.

    Raises
    ------
    KeyError
        If a named column is not in ``df``.
    ValueError
        If ``inner_radius`` is outside ``[0, 1)``, ``label`` is not one of the
        options, ``start`` is not finite, the palette is incomplete or a
        weight is negative.

    Notes
    -----
    Rows with a missing ``group_by`` or ``split_by`` are dropped with a warning
    that gives the count; without that a missing category would become its own
    grey slice. Rows with a missing weight are dropped by the stat with a
    warning. Slice angles are shares of each pie's own total, so a facet with
    fewer cells still closes. A pie with many levels is hard to read.

    Examples
    --------
    >>> import pandas as pd
    >>> import ggann as ag
    >>> df = pd.DataFrame({"kind": list("abcd"), "n": [4, 3, 2, 1]})
    >>> p = ag.ggpie(df, "kind", value="n", label="percent")
    """
    if label not in _LABELS:
        raise ValueError(f"label must be one of {_LABELS}, got {label!r}")
    if not 0 <= inner_radius < 1:
        raise ValueError(f"inner_radius must be in [0, 1), got {inner_radius!r}")
    if not np.isfinite(start):
        raise ValueError(f"start must be finite, got {start!r}")
    keep = [c for c in (group_by, value, split_by) if c is not None]
    frame = drop_missing_rows(df, [c for c in (group_by, split_by) if c is not None])
    frame, levels = prepare_levels(frame, group_by, order, keep=keep)

    mapping = aes(fill=group_by) if value is None else aes(fill=group_by, weight=value)
    plot = ggplot(frame, mapping) + geom_pie(color=border_color, size=border_size)
    if label != "none":
        text = after_stat("count") if label == "count" else after_stat("prop")
        extra = {"format_string": "{:g}" if label == "count" else "{:.0%}"}
        plot += geom_text(aes(label=text), stat=stat_pie, drop_zero=True, size=7, **extra)
    plot += coord_radial(theta="y", expand=False, inner_radius=inner_radius, start=start)
    if split_by is not None:
        plot += facet_wrap(f"~{split_by}")
    return (
        plot
        + palette_scales(palette, levels, ("fill",))
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


@accepts_scanpy_names
def ggring(
    df: pd.DataFrame,
    group_by: str,
    value: str | None = None,
    *,
    split_by: str | None = None,
    order=None,
    label: str = "none",
    inner_radius: float = 0.5,
    start: float = 0.0,
    palette=None,
    border_color: str = "white",
    border_size: float = 0.15,
):
    """Draw a ring (donut) of the shares of ``group_by`` in a data frame.

    Parameters
    ----------
    df : pandas.DataFrame
        One row per observation, or one row per category when ``value`` is
        given.
    group_by : str
        Column whose levels are the ring segments.
    value : str, optional
        Column of non-negative weights. ``None`` counts rows.
    split_by : str, optional
        Column that splits the data into one ring per level (facets).
    order : sequence, optional
        Segment order. ``None`` uses the categorical order of ``group_by``,
        then sorted values.
    label : {"none", "count", "percent"}
        Text drawn in each segment. Zero segments get no label.
    inner_radius : float
        Hole radius as a fraction of the outer radius, in ``(0, 1)``.
    start : float
        Angle of the first segment edge in radians, clockwise from 12 o'clock.
    palette : dict, list or str, optional
        A dict must cover every level. A list or a Matplotlib colormap name is
        used in level order. ``None`` keeps plotnine's default.
    border_color : str
        Colour of the segment outlines.
    border_size : float
        Width of the segment outlines.

    Returns
    -------
    plotnine.ggplot
        A plot that accepts further layers, scales and themes.

    Raises
    ------
    KeyError
        If a named column is not in ``df``.
    ValueError
        If ``inner_radius`` is outside ``(0, 1)``, ``label`` is not one of the
        options, ``start`` is not finite, the palette is incomplete or a
        weight is negative.

    Notes
    -----
    This is ``ggpie`` with a hole. Rows with a missing ``group_by`` or
    ``split_by`` are dropped with a warning that gives the count. See ``ggpie``
    for the other rules.

    Examples
    --------
    >>> import pandas as pd
    >>> import ggann as ag
    >>> df = pd.DataFrame({"kind": list("abcd"), "n": [4, 3, 2, 1]})
    >>> p = ag.ggring(df, "kind", value="n", inner_radius=0.6)
    """
    if not 0 < inner_radius < 1:
        raise ValueError(f"inner_radius must be in (0, 1) for a ring, got {inner_radius!r}")
    return ggpie(
        df,
        group_by,
        value,
        split_by=split_by,
        order=order,
        label=label,
        inner_radius=inner_radius,
        start=start,
        palette=palette,
        border_color=border_color,
        border_size=border_size,
    )
