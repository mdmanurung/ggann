"""Numpy geometry for the shape layers: Bezier curves, arcs and ribbons.

Private module with no plotnine import. Curves use the Bernstein form
``B(t) = sum_k C(m, k) t^k (1 - t)^(m - k) P_k`` with ``t_j = j / (n - 1)``.
Point arrays are ``float`` with shape ``(n, 2)`` unless a helper says it is
batched. Every helper with an ``n`` argument raises ``ValueError`` for ``n < 2``.
"""

from __future__ import annotations

from collections.abc import Sequence
from math import comb

import numpy as np
import pandas as pd

__all__ = [
    "arc_points",
    "bezier_points",
    "link_controls",
    "quad_control",
    "ribbon_polygon",
    "rows_from_pieces",
    "signed_area",
]


def _check_n(n: int) -> None:
    if n < 2:
        raise ValueError(f"n must be at least 2, got {n}")


def bezier_points(control: np.ndarray, n: int) -> np.ndarray:
    """Return ``n`` points along a Bezier curve of any degree.

    ``control`` is ``(m + 1, 2)`` or batched ``(F, m + 1, 2)``. The result is
    ``(n, 2)`` or ``(F, n, 2)``.
    """
    _check_n(n)
    control = np.asarray(control, dtype=float)
    m = control.shape[-2] - 1
    t = np.linspace(0.0, 1.0, n)[:, None]
    k = np.arange(m + 1)[None, :]
    weights = np.array([comb(m, j) for j in range(m + 1)], dtype=float)
    basis = weights * t**k * (1.0 - t) ** (m - k)
    return basis @ control


def quad_control(p0: np.ndarray, p1: np.ndarray, curvature: float) -> np.ndarray:
    """Return the control point of a quadratic curve from ``p0`` to ``p1``.

    ``C = (p0 + p1) / 2 - curvature * n_L`` where ``n_L = (-dy, dx)`` is the
    left normal of ``p1 - p0``. The apex is ``curvature / 2`` normals from the
    midpoint, so the bulge height is ``|curvature| * |p1 - p0| / 2``. Positive
    curvature bends right of travel; 0 is straight. ``p0`` and ``p1`` may be
    ``(2,)`` or batched ``(F, 2)``.
    """
    p0 = np.asarray(p0, dtype=float)
    p1 = np.asarray(p1, dtype=float)
    d = p1 - p0
    left = np.stack([-d[..., 1], d[..., 0]], axis=-1)
    return (p0 + p1) / 2.0 - curvature * left


def link_controls(
    p0: np.ndarray,
    p1: np.ndarray,
    strength: float = 0.5,
    axis: str = "x",
) -> np.ndarray:
    """Return the four control points of a cubic link with flat end tangents.

    For ``axis="x"`` the inner points are ``(x0 + s (x1 - x0), y0)`` and
    ``(x1 - s (x1 - x0), y1)``. ``axis="y"`` swaps the roles of x and y.
    ``p0`` and ``p1`` may be ``(2,)`` or batched ``(F, 2)``; the result is
    ``(4, 2)`` or ``(F, 4, 2)``.
    """
    if axis not in ("x", "y"):
        raise ValueError(f"axis must be 'x' or 'y', got {axis!r}")
    p0 = np.asarray(p0, dtype=float)
    p1 = np.asarray(p1, dtype=float)
    a = 0 if axis == "x" else 1
    step = strength * (p1[..., a] - p0[..., a])
    q1 = p0.copy()
    q2 = p1.copy()
    q1[..., a] = p0[..., a] + step
    q2[..., a] = p1[..., a] - step
    return np.stack([p0, q1, q2, p1], axis=-2)


def arc_points(
    theta0: float,
    theta1: float,
    r: float,
    n: int,
    center: tuple[float, float] = (0.0, 0.0),
) -> np.ndarray:
    """Return ``n`` points on a circular arc.

    Angles are radians, counter-clockwise from +x, linear in the point index.
    ``theta1 < theta0`` runs clockwise.
    """
    _check_n(n)
    theta = np.linspace(theta0, theta1, n)
    return np.column_stack([center[0] + r * np.cos(theta), center[1] + r * np.sin(theta)])


def ribbon_polygon(upper: np.ndarray, lower: np.ndarray) -> np.ndarray:
    """Join two edges into one open ring: ``upper`` forward, ``lower`` back."""
    return np.vstack([np.asarray(upper, dtype=float), np.asarray(lower, dtype=float)[::-1]])


def signed_area(ring: np.ndarray) -> float:
    """Return the shoelace area of a ring; positive is counter-clockwise."""
    ring = np.asarray(ring, dtype=float)
    nxt = np.roll(ring, -1, axis=0)
    return float(0.5 * np.sum(ring[:, 0] * nxt[:, 1] - nxt[:, 0] * ring[:, 1]))


def rows_from_pieces(pieces: Sequence[np.ndarray], keys: pd.DataFrame) -> pd.DataFrame:
    """Stack ``(n_i, 2)`` pieces into one long frame.

    ``keys`` has one row per piece. The result has ``x``, ``y``, ``group``
    (1..k) and the key columns repeated for each point of a piece. No ``t``
    column. An empty list gives an empty frame with the same columns.
    """
    if len(pieces) != len(keys):
        raise ValueError(f"got {len(pieces)} pieces for {len(keys)} key rows")
    if not pieces:
        out = pd.DataFrame({"x": [], "y": [], "group": np.array([], dtype=int)})
        return pd.concat([out, keys.iloc[:0].reset_index(drop=True)], axis=1)
    sizes = [len(p) for p in pieces]
    xy = np.vstack([np.asarray(p, dtype=float) for p in pieces])
    out = pd.DataFrame(
        {"x": xy[:, 0], "y": xy[:, 1], "group": np.repeat(np.arange(1, len(pieces) + 1), sizes)}
    )
    repeated = keys.reset_index(drop=True).loc[np.repeat(np.arange(len(pieces)), sizes)]
    return pd.concat([out, repeated.reset_index(drop=True)], axis=1)
