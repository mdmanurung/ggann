"""A hexagonal bin layer: ``stat_binhex``, ``geom_hex`` and ``gghex``."""

from __future__ import annotations

import math
from collections.abc import Callable, Sequence

import numpy as np
import pandas as pd
from plotnine import (
    aes,
    after_scale,
    after_stat,
    coord_equal,
    facet_wrap,
    geom_polygon,
    ggplot,
)
from plotnine.scales.scale_discrete import scale_discrete

from .._plotnine_compat import stat
from ._hexbin import _FUNS, _check_fun, bin_hexagons, hexagon_vertices

__all__ = ["geom_hex", "gghex", "stat_binhex"]

_STATS = ("count", "density", "ncount", "value")


def _pair(value, name: str, *, integer: bool) -> tuple:
    """``value`` as a validated pair; a scalar is a pair of one number."""
    items = tuple(value) if isinstance(value, (tuple, list, np.ndarray)) else (value,)
    if len(items) not in (1, 2):
        raise ValueError(f"{name} must be a number or a pair of numbers, got {value!r}")
    for item in items:
        if isinstance(item, bool) or not isinstance(item, (int, float, np.integer, np.floating)):
            raise ValueError(f"{name} must be numeric, got {value!r}")
        if integer and item != int(item):
            raise ValueError(f"{name} must be a whole number, got {value!r}")
        if not (math.isfinite(item) and item > 0):
            raise ValueError(f"{name} must be positive and finite, got {value!r}")
    return items


def _check_size(bins, binwidth) -> None:
    """Raise ``ValueError`` unless ``bins`` or ``binwidth`` is a valid size."""
    if bins is None and binwidth is None:
        raise ValueError(
            "give bins (a count across x) or binwidth (hexagon width in x units); "
            "there is no default"
        )
    if binwidth is not None:
        _pair(binwidth, "binwidth", integer=False)
    elif bins is not None:
        _pair(bins, "bins", integer=True)


def _cell_size(bins, binwidth, xdim, ydim) -> tuple[float, float]:
    """Hexagon ``(width, height)`` in scale units.

    ``binwidth`` wins over ``bins``. A scalar sets the width, and the height is
    ``2 * width / sqrt(3)``, which is regular under ``coord_equal``. A pair is
    ``(width, height)`` for ``binwidth`` and ``(nx, ny)`` for ``bins``, where
    width is ``x range / nx`` and height is ``y range / ny``.
    """
    if binwidth is not None:
        w, *rest = _pair(binwidth, "binwidth", integer=False)
        return float(w), float(rest[0]) if rest else 2 * float(w) / math.sqrt(3)
    nx, *rest = _pair(bins, "bins", integer=True)
    w = (xdim[1] - xdim[0]) / nx
    if rest:
        return w, (ydim[1] - ydim[0]) / rest[0]
    return w, 2 * w / math.sqrt(3)


class stat_binhex(stat):
    """Count points in hexagons, then emit six polygon vertices per bin.

    Hexagons are pointy-top, found by nearest centre on a triangular lattice.
    The lattice origin is the lower end of the x and y scale limits, so panels
    with fixed scales share one lattice. Binning runs after scale transforms,
    so a log axis bins in log units and ``binwidth`` is in those units.

    Parameters
    ----------
    mapping : plotnine.mapping.aes, optional
        Aesthetics. ``x`` and ``y`` are required. ``weight`` sums weights
        instead of counting points. ``z`` is summarised per bin by ``fun``.
    data : pandas.DataFrame, optional
        Layer data. Defaults to the plot data.
    **kwargs
        ``bins`` is an int (hexagons across x; the hexagon is regular under
        ``coord_equal``) or a pair ``(nx, ny)`` (width is the x range over
        ``nx``, height the y range over ``ny``). ``binwidth`` is a number
        (width in x units, height ``2 * width / sqrt(3)``) or a pair
        ``(width, height)``, and wins over ``bins``. One of the two is
        required; there is no default. ``fun`` is ``"mean"``, ``"median"``,
        ``"sum"``, ``"min"``, ``"max"`` or a callable; it is required when
        ``z`` is mapped, and has no default. ``geom`` is the geom that draws
        the polygons (default ``geom_polygon``). Other keywords are constant
        aesthetics.

    Returns
    -------
    stat_binhex
        A layer to add to a ``plotnine.ggplot``. Its output has six rows per
        bin with ``x``, ``y`` (the vertex), ``group`` and ``hex`` (bin id 1 to
        k in the panel), ``i``, ``j`` (lattice index), ``count`` (sum of
        weights), ``density`` (``count / sum(count)``), ``ncount`` (``count /
        max(count)``), ``value`` (only when ``z`` is mapped) and ``PANEL``.

    Raises
    ------
    ValueError
        When the plot is built, if neither ``bins`` nor ``binwidth`` is given,
        either is not positive, ``z`` is mapped without a valid ``fun``, ``z``
        or ``weight`` holds a missing value, or the resolved hexagon size is
        not positive and finite.
    TypeError
        If x or y uses a discrete scale.
    plotnine.exceptions.PlotnineError
        If a keyword matches no stat, geom or aesthetic parameter.

    Notes
    -----
    ``group`` is set by the stat and replaces any mapped ``group``; the stat
    pools every input group in a panel, so a mapped discrete ``color`` is
    dropped. Rows with a missing or non-finite ``x`` or ``y`` are removed with
    a warning (``na_rm=True`` silences it), and scale limits turn out-of-range
    points into missing ones. ``count`` always counts the points; a missing
    ``z`` is an error, not skipped. With free facet scales each panel has its
    own lattice.

    Examples
    --------
    >>> import numpy as np
    >>> import pandas as pd
    >>> from plotnine import aes, ggplot
    >>> import ggann as ag
    >>> df = pd.DataFrame({"x": np.arange(50.0), "y": np.arange(50.0) ** 0.5})
    >>> p = ggplot(df, aes("x", "y")) + ag.stat_binhex(bins=10)
    """

    REQUIRED_AES = {"x", "y"}
    DEFAULT_AES = {"fill": after_stat("count"), "weight": None, "z": None}
    DEFAULT_PARAMS = {"geom": geom_polygon, "bins": None, "binwidth": None, "fun": None}
    CREATES = {"count", "density", "ncount", "value", "hex"}

    def setup_params(self, data):
        _check_size(self.params["bins"], self.params["binwidth"])
        if "z" in data.columns:
            if self.params["fun"] is None:
                raise ValueError(
                    f"z is mapped but fun is not set; use one of {list(_FUNS)} or a callable"
                )
            _check_fun(self.params["fun"])
        return self.params

    def compute_panel(self, data, scales):
        if isinstance(scales.x, scale_discrete) or isinstance(scales.y, scale_discrete):
            raise TypeError("stat_binhex needs continuous x and y scales, not discrete ones")
        has_z = "z" in data.columns
        panel = data["PANEL"].iloc[0] if len(data) else None
        columns = ["x", "y", "group", "hex", "i", "j", "count", "density", "ncount"]
        if has_z:
            columns.append("value")
        if len(data) == 0:
            out = pd.DataFrame({c: pd.Series(dtype=float) for c in columns})
            out[["group", "hex", "i", "j"]] = out[["group", "hex", "i", "j"]].astype(np.int64)
            out["PANEL"] = pd.Series([], dtype=data["PANEL"].dtype)
            return out

        xdim, ydim = scales.x.dimension(), scales.y.dimension()
        w, h = _cell_size(self.params["bins"], self.params["binwidth"], xdim, ydim)
        if not (np.isfinite(w) and np.isfinite(h) and w > 0 and h > 0):
            raise ValueError(f"resolved hexagon size must be positive and finite, got {w}, {h}")
        cells = bin_hexagons(
            data["x"],
            data["y"],
            width=w,
            height=h,
            origin=(xdim[0], ydim[0]),
            weight=data["weight"] if "weight" in data.columns else None,
            z=data["z"] if has_z else None,
            fun=self.params["fun"],
        )
        k = len(cells)
        vx, vy = hexagon_vertices(cells["x"], cells["y"], w, h)
        ids = np.arange(1, k + 1)
        count = cells["count"].to_numpy()
        out = pd.DataFrame(
            {
                "x": vx.ravel(),
                "y": vy.ravel(),
                "group": np.repeat(ids, 6),
                "hex": np.repeat(ids, 6),
                "i": np.repeat(cells["i"].to_numpy(), 6),
                "j": np.repeat(cells["j"].to_numpy(), 6),
                "count": np.repeat(count, 6),
                "density": np.repeat(count / count.sum(), 6),
                "ncount": np.repeat(count / count.max(), 6),
            }
        )
        if has_z:
            out["value"] = np.repeat(cells["value"].to_numpy(), 6)
        out["PANEL"] = pd.Series([panel] * len(out), dtype=data["PANEL"].dtype)
        return out


class geom_hex(geom_polygon):
    """Draw hexagonal bins of ``x`` and ``y``, filled by their count.

    A polygon geom with ``stat_binhex`` as its default stat. Unless you give a
    ``color``, the outline copies the fill at width 0.2, which hides the white
    hairlines between hexagons. That rule is skipped when ``alpha`` is below 1
    (an opaque outline would draw a grid) or mapped, and a constant ``fill``
    becomes the outline colour.

    Parameters
    ----------
    mapping : plotnine.mapping.aes, optional
        Aesthetics. ``x`` and ``y`` are required; ``weight`` and ``z`` are
        optional; ``fill`` defaults to ``after_stat("count")``.
    data : pandas.DataFrame, optional
        Layer data. Defaults to the plot data.
    **kwargs
        ``bins`` or ``binwidth`` (one is required) and ``fun`` (required with
        ``z``) go to ``stat_binhex``. ``raster`` and ``alpha`` go to the
        polygon geom. Other keywords are constant aesthetics.

    Returns
    -------
    geom_hex
        A layer to add to a ``plotnine.ggplot``.

    Raises
    ------
    ValueError
        When the plot is built, if neither ``bins`` nor ``binwidth`` is given
        or a size is not positive; see ``stat_binhex`` for the rest.
    TypeError
        If x or y uses a discrete scale.
    plotnine.exceptions.PlotnineError
        If a keyword matches no stat, geom or aesthetic parameter.

    Notes
    -----
    One polygon is drawn per bin, so vector output grows with the bin count;
    ``raster=True`` rasterises the layer. Use ``coord_equal()`` for regular
    hexagons with a scalar ``bins`` or ``binwidth``. See ``stat_binhex`` for
    how bins are placed and for the missing-value rules.

    Examples
    --------
    >>> import numpy as np
    >>> import pandas as pd
    >>> from plotnine import aes, ggplot, coord_equal
    >>> import ggann as ag
    >>> df = pd.DataFrame({"x": np.arange(50.0), "y": np.arange(50.0) ** 0.5})
    >>> p = ggplot(df, aes("x", "y")) + ag.geom_hex(bins=10) + coord_equal()
    """

    DEFAULT_PARAMS = {"stat": stat_binhex}
    DEFAULT_AES = {**geom_polygon.DEFAULT_AES, "size": 0.2}

    def __init__(self, mapping=None, data=None, **kwargs):
        mapped = set(mapping) if mapping is not None else set()
        given = set(kwargs) | mapped
        opaque = "alpha" not in mapped and kwargs.get("alpha", 1) == 1
        if opaque and not ({"color", "colour"} & given):
            if "fill" in kwargs:
                kwargs["color"] = kwargs["fill"]
            else:
                mapping = aes(**{**dict(mapping or {}), "color": after_scale("fill")})
        super().__init__(mapping, data, **kwargs)


def gghex(
    df: pd.DataFrame,
    x: str,
    y: str,
    *,
    z: str | None = None,
    fun: str | Callable | None = None,
    weight: str | None = None,
    bins: int | Sequence[int] | None = None,
    binwidth: float | Sequence[float] | None = None,
    stat: str = "count",
    facet: str | None = None,
    equal_aspect: bool = True,
):
    """Hexagonal bin plot of two columns of a data frame.

    Parameters
    ----------
    df : pandas.DataFrame
        Long-format table.
    x, y : str
        Numeric columns.
    z : str, optional
        Numeric column summarised per hexagon by ``fun``.
    fun : {"mean", "median", "sum", "min", "max"} or callable, optional
        Summary of ``z`` per hexagon. Required with ``z``; there is no default.
    weight : str, optional
        Numeric column; the count is the sum of weights.
    bins : int or pair of int, optional
        Hexagons across x, or ``(nx, ny)``. See :class:`ggann.stat_binhex`.
    binwidth : float or pair of float, optional
        Hexagon width in scale units, or ``(width, height)``. Wins over ``bins``.
        One of ``bins`` and ``binwidth`` is required; there is no default.
    stat : {"count", "density", "ncount", "value"}, default="count"
        Bin statistic that sets the fill. ``"value"`` needs ``z``.
    facet : str, optional
        Column to facet by (``facet_wrap``). Panels share one lattice.
    equal_aspect : bool, default=True
        Add ``coord_equal()``, so a scalar ``bins`` gives regular hexagons.

    Returns
    -------
    plotnine.ggplot
        A composable plot with one ``geom_hex`` layer.

    Raises
    ------
    KeyError
        If a named column is missing from ``df``.
    ValueError
        For an unknown ``stat`` or ``fun``, ``stat="value"`` without ``z``,
        ``z`` without ``fun``, neither ``bins`` nor ``binwidth``, or a size
        that is not positive.

    Notes
    -----
    Rows with a missing ``x`` or ``y`` are dropped with a warning when the plot
    is built. A missing ``z`` or ``weight`` raises ``ValueError`` then. The
    table stays per-row; ``plot.layer_data(0)`` holds the hexagons.

    Examples
    --------
    >>> import numpy as np
    >>> import pandas as pd
    >>> import ggann as ag
    >>> rng = np.random.default_rng(0)
    >>> df = pd.DataFrame({"x": rng.normal(size=500), "y": rng.normal(size=500)})
    >>> p = ag.gghex(df, "x", "y", bins=20)
    """
    from ..publication import _family_theme

    if stat not in _STATS:
        raise ValueError(f"stat must be one of {list(_STATS)}, got {stat!r}")
    if stat == "value" and z is None:
        raise ValueError("stat='value' needs z")
    if z is not None:
        if fun is None:
            raise ValueError(f"z needs fun: one of {list(_FUNS)} or a callable")
        _check_fun(fun)
    _check_size(bins, binwidth)
    named = {"x": x, "y": y, "z": z, "weight": weight, "facet": facet}
    for role, column in named.items():
        if column is not None and column not in df.columns:
            raise KeyError(f"{role} column {column!r} is not in df")

    mapping = aes(x, y, fill=after_stat(stat))
    if z is not None:
        mapping["z"] = z
    if weight is not None:
        mapping["weight"] = weight
    layer = geom_hex(mapping, bins=bins, binwidth=binwidth, fun=fun)
    plot = ggplot(df) + layer + _family_theme("standard")
    if equal_aspect:
        plot += coord_equal()
    if facet is not None:
        plot += facet_wrap("~" + facet)
    return plot
