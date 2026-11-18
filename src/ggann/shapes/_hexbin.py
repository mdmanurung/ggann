"""Hexagonal binning in numpy: nearest centre on a triangular lattice.

Private module with no plotnine import. Hexagons are pointy-top. ``width`` is
the flat-to-flat width in x units and ``height`` the vertex-to-vertex height in
y units. Centre rows sit ``0.75 * height`` apart and odd rows shift by
``width / 2``. A point belongs to the hexagon whose centre is nearest in the
metric where a regular hexagon is a Voronoi cell (Carr et al. 1987, JASA 82:424).
"""

from __future__ import annotations

from collections.abc import Callable, Sequence

import numpy as np
import pandas as pd

__all__ = ["bin_hexagons", "hexagon_vertices"]

_FUNS = ("mean", "median", "sum", "min", "max")


def _check_fun(fun) -> None:
    if callable(fun) or fun in _FUNS:
        return
    raise ValueError(f"fun must be one of {list(_FUNS)} or a callable, got {fun!r}")


def _assign(x: np.ndarray, y: np.ndarray, width: float, height: float, origin):
    """Lattice index ``(j, i)`` of the nearest hexagon centre for each point."""
    x0, y0 = origin
    u = (x - x0) / width
    v = (y - y0) / (0.75 * height)
    j0 = np.floor(v).astype(np.int64)
    best_j = j0
    best_i = np.floor(u - 0.5 * (j0 % 2) + 0.5).astype(np.int64)
    best_d = (u - best_i - 0.5 * (j0 % 2)) ** 2 + 0.75 * (v - j0) ** 2
    # Two candidate rows suffice: a row's hexagons span v in [j - 2/3, j + 2/3].
    j1 = j0 + 1
    off = 0.5 * (j1 % 2)
    i1 = np.floor(u - off + 0.5).astype(np.int64)
    d1 = (u - i1 - off) ** 2 + 0.75 * (v - j1) ** 2
    take = d1 < best_d  # a tie keeps row j0
    return np.where(take, j1, best_j), np.where(take, i1, best_i)


def bin_hexagons(
    x: Sequence[float],
    y: Sequence[float],
    *,
    width: float,
    height: float,
    origin: tuple[float, float],
    weight: Sequence[float] | None = None,
    z: Sequence[float] | None = None,
    fun: str | Callable | None = None,
) -> pd.DataFrame:
    """Count points per hexagon.

    Returns one row per non-empty bin, sorted by ``(j, i)``, with columns ``j``,
    ``i``, ``x``, ``y`` (the centre) and ``count`` (the sum of ``weight``, or
    the number of points), plus ``value`` (``fun`` of ``z`` per bin) when ``z``
    is given. ``fun`` has no default: it is required with ``z``. Raises
    ``ValueError`` for a non-finite ``x``, ``y``, ``weight`` or ``z``, a
    non-positive ``width`` or ``height``, or ``z`` without a valid ``fun``.
    """
    xs = np.asarray(x, dtype=float)
    ys = np.asarray(y, dtype=float)
    if not (np.isfinite(width) and np.isfinite(height) and width > 0 and height > 0):
        raise ValueError(f"width and height must be positive and finite, got {width}, {height}")
    if not (np.isfinite(xs).all() and np.isfinite(ys).all()):
        raise ValueError("x and y must be finite; drop missing values before binning")
    w = np.ones(len(xs)) if weight is None else np.asarray(weight, dtype=float)
    if not np.isfinite(w).all():
        raise ValueError("weight must be finite")
    if z is not None:
        if fun is None:
            raise ValueError(f"z needs fun: one of {list(_FUNS)} or a callable")
        _check_fun(fun)
        zv = np.asarray(z, dtype=float)
        if not np.isfinite(zv).all():
            raise ValueError("z must be finite; drop or fill missing values first")
    columns = ["j", "i", "x", "y", "count"] + (["value"] if z is not None else [])
    if len(xs) == 0:
        return pd.DataFrame({c: pd.Series(dtype=float) for c in columns}).astype(
            {"j": np.int64, "i": np.int64}
        )

    j, i = _assign(xs, ys, width, height, origin)
    frame = pd.DataFrame({"j": j, "i": i, "count": w})
    if z is not None:
        frame["z"] = zv
    grouped = frame.groupby(["j", "i"], sort=True)
    out = grouped["count"].sum().reset_index()
    if z is not None:
        out["value"] = grouped["z"].agg(fun).to_numpy()
    out["x"] = origin[0] + width * (out["i"] + 0.5 * (out["j"] % 2))
    out["y"] = origin[1] + 0.75 * height * out["j"]
    return pd.DataFrame(out[columns])


def hexagon_vertices(
    cx: np.ndarray, cy: np.ndarray, width: float, height: float
) -> tuple[np.ndarray, np.ndarray]:
    """Six vertices per centre, clockwise from the top; arrays of shape ``(k, 6)``."""
    cx = np.asarray(cx, dtype=float)[:, None]
    cy = np.asarray(cy, dtype=float)[:, None]
    dx = np.array([0, 0.5, 0.5, 0, -0.5, -0.5]) * width
    dy = np.array([0.5, 0.25, -0.25, -0.5, -0.25, 0.25]) * height
    return cx + dx, cy + dy
