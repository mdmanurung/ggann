"""A Bezier curve layer: one bent line per row of ``x, y, xend, yend``."""

from __future__ import annotations

import numpy as np
import pandas as pd
from plotnine import geom_path

from .._plotnine_compat import stat
from ._geometry import bezier_points, quad_control

__all__ = ["geom_curve", "stat_curve"]


class stat_curve(stat):
    """Replace each ``(x, y) -> (xend, yend)`` row by ``n`` points on a curve.

    The curve is a quadratic Bezier. Its apex sits ``curvature / 2`` normals
    from the segment midpoint, so the bulge is ``|curvature| / 2`` of the
    chord length. Positive ``curvature`` bends right of the direction of
    travel, negative bends left, 0 is a straight line.

    Parameters
    ----------
    mapping : plotnine.mapping.aes, optional
        Aesthetics. ``x``, ``y``, ``xend`` and ``yend`` are required.
    data : pandas.DataFrame, optional
        Layer data. Defaults to the plot data.
    **kwargs
        ``curvature`` (float, default 0.5) sets the bend. ``n`` (int, default
        50) is the number of points per curve, at least 2. ``geom`` is the
        geom that draws the points (default ``geom_path``). ``na_rm`` and
        ``position`` work as in any plotnine stat. Other keywords are
        constant aesthetics.

    Returns
    -------
    stat_curve
        A layer to add to a ``plotnine.ggplot``. Its output has ``x``, ``y``,
        ``group``, ``t`` in ``[0, 1]``, ``PANEL`` and every other input
        column, ordered by input row and then ``t``. ``xend`` and ``yend``
        are dropped.

    Raises
    ------
    ValueError
        When the plot is built, if ``n`` is not an integer of at least 2 or
        ``curvature`` is not finite.
    plotnine.exceptions.PlotnineError
        If a keyword matches no stat, geom or aesthetic parameter.

    Notes
    -----
    ``group`` is set by the stat: it numbers the rows that survive the
    missing-value filter, 1 to k within each panel, and replaces any mapped
    ``group``. Rows with a missing or non-finite endpoint are removed with a
    warning (``na_rm=True`` silences it).

    Curvature is measured in transformed scale units, so the shape is exact
    only when x and y share units, as under ``coord_equal``. With unequal
    units the bulge follows the larger span: a segment from ``(0, 0)`` to
    ``(1, 1000)`` at the default curvature reaches x = 250.5. Curves are not
    supported under ``coord_radial``, where the path bends in (theta, r)
    space and can push r below 0. Under ``coord_flip`` the bend is assumed to
    mirror.

    Once ggann is imported, ``stat="curve"`` in any plotnine layer resolves to
    this class, because plotnine keeps one registry keyed by class name.

    Examples
    --------
    >>> import pandas as pd
    >>> from plotnine import aes, ggplot, geom_point
    >>> import ggann as ag
    >>> df = pd.DataFrame({"x": [0, 0], "y": [0, 0], "xend": [4, 0], "yend": [0, 2]})
    >>> p = ggplot(df, aes("x", "y", xend="xend", yend="yend")) + ag.stat_curve(n=20)
    """

    REQUIRED_AES = {"x", "y", "xend", "yend"}
    DEFAULT_PARAMS = {"geom": geom_path, "curvature": 0.5, "n": 50}
    DROPPED_AES = ["xend", "yend"]
    CREATES = {"t"}

    def setup_params(self, data):
        n = self.params["n"]
        curvature = self.params["curvature"]
        if isinstance(n, bool) or not isinstance(n, (int, np.integer)) or n < 2:
            raise ValueError(f"n must be an integer of at least 2, got {n!r}")
        if not np.isfinite(curvature):
            raise ValueError(f"curvature must be finite, got {curvature!r}")
        return self.params

    def compute_panel(self, data, scales):
        n = int(self.params["n"])
        rest = data.drop(columns=["xend", "yend"])
        if len(data) == 0:
            return rest.assign(t=pd.Series(dtype=float), group=pd.Series(dtype=int))

        p0 = data[["x", "y"]].to_numpy(dtype=float)
        p1 = data[["xend", "yend"]].to_numpy(dtype=float)
        control = np.stack([p0, quad_control(p0, p1, self.params["curvature"]), p1], axis=1)
        xy = bezier_points(control, n).reshape(-1, 2)

        out = rest.iloc[np.repeat(np.arange(len(data)), n)].reset_index(drop=True)
        out["x"] = xy[:, 0]
        out["y"] = xy[:, 1]
        out["t"] = np.tile(np.linspace(0.0, 1.0, n), len(data))
        out["group"] = np.repeat(np.arange(1, len(data) + 1), n)
        return out


class geom_curve(geom_path):
    """Draw a quadratic Bezier curve for each ``(x, y) -> (xend, yend)`` row.

    A path geom with ``stat_curve`` as its default stat. It takes the path
    aesthetics (``alpha``, ``color``, ``linetype``, ``size``) and the path
    parameters ``lineend``, ``linejoin`` and ``arrow``.

    Parameters
    ----------
    mapping : plotnine.mapping.aes, optional
        Aesthetics. ``x``, ``y``, ``xend`` and ``yend`` are required.
    data : pandas.DataFrame, optional
        Layer data. Defaults to the plot data.
    **kwargs
        ``curvature`` (float, default 0.5) and ``n`` (int, default 50) go to
        ``stat_curve``. ``arrow``, ``lineend`` and ``linejoin`` go to the
        path geom. Other keywords are constant aesthetics.

    Returns
    -------
    geom_curve
        A layer to add to a ``plotnine.ggplot``.

    Raises
    ------
    ValueError
        When the plot is built, if ``n`` is not an integer of at least 2 or
        ``curvature`` is not finite.
    plotnine.exceptions.PlotnineError
        If a keyword matches no stat, geom or aesthetic parameter.

    Notes
    -----
    Positive ``curvature`` bends right of the direction of travel. Curvature
    is measured in transformed scale units, so the shape is exact only when x
    and y share units, as under ``coord_equal``. A mapped ``group`` is
    ignored, rows with missing endpoints are removed with a warning, and
    curves are not supported under ``coord_radial``. See ``stat_curve`` for
    the details.

    Examples
    --------
    >>> import pandas as pd
    >>> from plotnine import aes, ggplot
    >>> import ggann as ag
    >>> df = pd.DataFrame({"x": [0, 0], "y": [0, 0], "xend": [4, 0], "yend": [0, 2]})
    >>> p = ggplot(df, aes("x", "y", xend="xend", yend="yend")) + ag.geom_curve(curvature=0.5)
    """

    DEFAULT_PARAMS = {"stat": stat_curve}
