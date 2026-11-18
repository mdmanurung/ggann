"""Hand-recomputed checks for the private shape geometry helpers."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from ggann.shapes._geometry import (
    arc_points,
    bezier_points,
    link_controls,
    quad_control,
    ribbon_polygon,
    rows_from_pieces,
    signed_area,
)


def test_quad_control_bends_right_of_travel():
    assert quad_control((0, 0), (4, 0), 0.5) == pytest.approx([2, -2])
    assert quad_control((0, 0), (0, 2), 1) == pytest.approx([2, 1])


def test_quad_control_zero_curvature_is_the_midpoint():
    assert quad_control((0, 0), (4, 2), 0) == pytest.approx([2, 1])


def test_quadratic_bezier_points():
    out = bezier_points([(0, 0), (2, -2), (4, 0)], 5)
    assert out[:, 0] == pytest.approx([0, 1, 2, 3, 4])
    assert out[:, 1] == pytest.approx([0, -0.75, -1, -0.75, 0])


def test_zero_curvature_stays_on_the_segment():
    control = np.stack([[0, 0], quad_control((0, 0), (4, 2), 0), [4, 2]])
    out = bezier_points(control, 7)
    assert out[:, 1] == pytest.approx(out[:, 0] / 2)


def test_cubic_bezier_points_and_link_controls():
    control = link_controls((0, 0), (1, 1))
    assert control == pytest.approx(np.array([(0, 0), (0.5, 0), (0.5, 1), (1, 1)]))
    out = bezier_points(control, 5)
    expected = [(0, 0), (0.296875, 0.15625), (0.5, 0.5), (0.703125, 0.84375), (1, 1)]
    assert out == pytest.approx(np.array(expected))


def test_link_controls_axis_y_swaps_roles():
    control = link_controls((0, 0), (1, 2), axis="y")
    assert control == pytest.approx(np.array([(0, 0), (0, 1), (1, 1), (1, 2)]))
    with pytest.raises(ValueError, match="axis"):
        link_controls((0, 0), (1, 1), axis="z")


def test_arc_points_quarter_turn():
    out = arc_points(0, np.pi / 2, 1, 3)
    assert out == pytest.approx(np.array([(1, 0), (0.70710678, 0.70710678), (0, 1)]))


def test_arc_points_clockwise_when_theta1_is_smaller():
    out = arc_points(np.pi / 2, 0, 2, 3, center=(1, 1))
    assert out[0] == pytest.approx([1, 3])
    assert out[-1] == pytest.approx([3, 1])


def test_ribbon_polygon_and_signed_area():
    upper = np.array([(0, 1), (0.5, 1.5), (1, 2)])
    lower = np.array([(0, 0), (0.5, 0.5), (1, 1)])
    ring = ribbon_polygon(upper, lower)
    assert len(ring) == 6
    assert ring[-1] == pytest.approx([0, 0])
    assert signed_area(ring) == pytest.approx(-1.0)
    assert signed_area(np.array([(0, 0), (1, 0), (1, 1), (0, 1)])) == pytest.approx(1.0)


def test_point_count_below_two_raises():
    with pytest.raises(ValueError, match="at least 2"):
        bezier_points([(0, 0), (1, 1)], 1)
    with pytest.raises(ValueError, match="at least 2"):
        arc_points(0, 1, 1, 1)


def test_batched_calls_equal_single_calls():
    control = np.array([[(0, 0), (2, -2), (4, 0)], [(0, 0), (1, 3), (0, 2)]], dtype=float)
    batched = bezier_points(control, 6)
    assert batched.shape == (2, 6, 2)
    for i in range(2):
        assert batched[i] == pytest.approx(bezier_points(control[i], 6))

    p0 = np.array([(0, 0), (1, 1)], dtype=float)
    p1 = np.array([(4, 0), (1, 3)], dtype=float)
    quad = quad_control(p0, p1, 0.5)
    link = link_controls(p0, p1)
    assert link.shape == (2, 4, 2)
    for i in range(2):
        assert quad[i] == pytest.approx(quad_control(p0[i], p1[i], 0.5))
        assert link[i] == pytest.approx(link_controls(p0[i], p1[i]))


def test_rows_from_pieces_repeats_the_keys():
    keys = pd.DataFrame({"kind": ["a", "b"]})
    out = rows_from_pieces([np.zeros((2, 2)), np.ones((3, 2))], keys)
    assert list(out.columns) == ["x", "y", "group", "kind"]
    assert out["group"].tolist() == [1, 1, 2, 2, 2]
    assert out["kind"].tolist() == ["a", "a", "b", "b", "b"]
    assert out["x"].tolist() == [0, 0, 1, 1, 1]


def test_rows_from_pieces_empty_keeps_columns():
    keys = pd.DataFrame({"kind": ["a"]})
    out = rows_from_pieces([], keys.iloc[:0])
    assert len(out) == 0
    assert {"x", "y", "group", "kind"} <= set(out.columns)
    with pytest.raises(ValueError, match="key rows"):
        rows_from_pieces([np.zeros((2, 2))], keys.iloc[:0])
