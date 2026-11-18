"""stat_curve and geom_curve: tables, edge cases, errors and renders."""

from __future__ import annotations

import warnings
from io import BytesIO

import numpy as np
import pandas as pd
import plotnine as p9
import pytest
from plotnine.exceptions import PlotnineError

import ggann as ag
import ggann.shapes as shapes


def _frame() -> pd.DataFrame:
    return pd.DataFrame({"x": [0, 0], "y": [0, 0], "xend": [4, 0], "yend": [0, 2], "k": ["a", "b"]})


def _plot(df, layer, **aes_extra) -> p9.ggplot:
    return p9.ggplot(df, p9.aes("x", "y", xend="xend", yend="yend", **aes_extra)) + layer


def _render(plot) -> None:
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        plot.save(BytesIO(), format="png", verbose=False)


def test_exported_from_top_level_and_shapes():
    assert ag.geom_curve is shapes.geom_curve
    assert ag.stat_curve is shapes.stat_curve
    assert {"geom_curve", "stat_curve"} <= set(ag.__all__)


def test_stat_table_for_two_rows():
    out = _plot(_frame(), ag.geom_curve(p9.aes(color="k"), n=5)).layer_data(0)
    assert len(out) == 10
    assert out["group"].tolist() == [1] * 5 + [2] * 5
    assert out["t"].tolist() == pytest.approx([0, 0.25, 0.5, 0.75, 1] * 2)
    g2 = out[out["group"] == 2]
    assert g2["x"].tolist() == pytest.approx([0, 0.375, 0.5, 0.375, 0])
    assert g2["y"].tolist() == pytest.approx([0, 0.5, 1, 1.5, 2])
    assert "xend" not in out.columns
    assert out.groupby("group")["color"].nunique().tolist() == [1, 1]
    assert out["color"].nunique() == 2


def test_negative_curvature_bends_left():
    out = _plot(_frame(), ag.geom_curve(curvature=-0.3, n=5)).layer_data(0)
    assert out.loc[out["group"] == 1, "y"].tolist() == pytest.approx([0, 0.45, 0.6, 0.45, 0])
    assert out.loc[out["group"] == 2, "x"].tolist() == pytest.approx([0, -0.225, -0.3, -0.225, 0])


def test_stat_layer_curvature_one():
    out = _plot(_frame(), ag.stat_curve(curvature=1, n=5)).layer_data(0)
    assert out.loc[out["group"] == 1, "y"].tolist() == pytest.approx([0, -1.5, -2, -1.5, 0])


def test_facets_restart_group_numbers():
    df = pd.DataFrame(
        {
            "x": [0, 0, 1],
            "y": [0, 0, 1],
            "xend": [1, 2, 3],
            "yend": [1, 2, 0],
            "f": ["u", "u", "v"],
        }
    )
    plot = _plot(df, ag.geom_curve(n=3)) + p9.facet_wrap("f")
    out = plot.layer_data(0)
    assert out.groupby("PANEL", observed=True)["group"].unique().map(list).tolist() == [
        [1, 2],
        [1],
    ]
    empty = _plot(df, ag.geom_curve(n=3)) + p9.facet_wrap("f", drop=False)
    _render(empty)


def test_missing_endpoints_are_dropped_with_one_warning():
    df = pd.DataFrame(
        {
            "x": [0, 0, 0, 0],
            "y": [0, 0, 0, 0],
            "xend": [1, np.nan, 1, np.inf],
            "yend": [1, 1, 2, 1],
        }
    )
    with pytest.warns(UserWarning, match="Removed 2 rows") as record:
        out = _plot(df, ag.geom_curve(n=3)).layer_data(0)
    assert len(record) == 1
    assert len(out) == 6
    assert sorted(out["group"].unique()) == [1, 2]
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        quiet = _plot(df, ag.geom_curve(n=3, na_rm=True)).layer_data(0)
    assert quiet["x"].tolist() == pytest.approx(out["x"].tolist())


def test_mapped_group_is_overwritten():
    out = _plot(_frame(), ag.geom_curve(n=3), group="k").layer_data(0)
    assert out["group"].tolist() == [1] * 3 + [2] * 3
    df = _frame().assign(k=["same", "same"])
    out = _plot(df, ag.geom_curve(n=3), group="k").layer_data(0)
    assert out["group"].tolist() == [1] * 3 + [2] * 3


def test_empty_frame_and_zero_length_segment_render_quietly():
    empty = _frame().iloc[:0]
    _render(_plot(empty, ag.geom_curve(n=3)))
    zero = pd.DataFrame({"x": [1.0], "y": [1.0], "xend": [1.0], "yend": [1.0]})
    _render(_plot(zero, ag.geom_curve(n=3, arrow=p9.arrow())))


def test_unequal_axis_units_stretch_the_x_range():
    df = pd.DataFrame({"x": [0.0], "y": [0.0], "xend": [1.0], "yend": [1000.0]})
    out = _plot(df, ag.geom_curve(n=5)).layer_data(0)
    assert out["x"].max() == pytest.approx(250.5)


def test_reversed_y_scale_keeps_the_bend_on_screen():
    out = (_plot(_frame(), ag.geom_curve(n=5)) + p9.scale_y_reverse()).layer_data(0)
    assert out.loc[out["group"] == 1, "y"].tolist() == pytest.approx([0, -0.75, -1, -0.75, 0])


def test_registry_string_resolves_to_the_stat():
    plot = _plot(_frame(), p9.geom_path(stat="curve", n=5))
    fig = plot.draw()
    assert len(fig.axes[0].lines) == 2
    assert plot.layer_data(0)["group"].nunique() == 2


@pytest.mark.parametrize("n", [1, 0, 2.5, True])
def test_bad_n_raises_at_build(n):
    with pytest.raises(ValueError, match="n must be"):
        _plot(_frame(), ag.geom_curve(n=n)).layer_data(0)


@pytest.mark.parametrize("curvature", [float("nan"), float("inf")])
def test_non_finite_curvature_raises_at_build(curvature):
    with pytest.raises(ValueError, match="curvature must be finite"):
        _plot(_frame(), ag.geom_curve(curvature=curvature)).layer_data(0)


def test_misspelt_keyword_raises():
    with pytest.raises(PlotnineError):
        _plot(_frame(), ag.geom_curve(curvatur=1))


def test_arrow_and_equal_coord_render_quietly():
    _render(_plot(_frame(), ag.geom_curve(arrow=p9.arrow())))
    _render(_plot(_frame(), ag.geom_curve()) + p9.coord_equal())


def test_drawn_arrow_params_reach_the_geom():
    fig = _plot(_frame(), ag.geom_curve(curvature=-0.3, n=5, arrow=p9.arrow())).draw()
    ax = fig.axes[0]
    assert len(ax.lines) == 2
    assert all(len(line.get_xdata()) == 5 for line in ax.lines)
