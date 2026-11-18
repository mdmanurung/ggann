"""bin_hexagons, stat_binhex, geom_hex, gghex and plot_embedding(hex=True)."""

from __future__ import annotations

import warnings
from io import BytesIO

import numpy as np
import pandas as pd
import plotnine as p9
import pytest
from plotnine.exceptions import PlotnineError, PlotnineWarning

import ggann as ag
import ggann.shapes as shapes
from ggann.shapes._hexbin import bin_hexagons, hexagon_vertices

X = [0.1, 0.2, 0.9, 1.1, 0.5, 0.45, 1.6, 2.0, 0.0, 0.5]
Y = [0.1, 0.0, 0.05, 0.0, 0.75, 0.80, 0.70, 1.5, 0.375, 0.375]
CELLS = [(0, 0), (0, 0), (1, 0), (1, 0), (0, 1), (0, 1), (1, 1), (2, 2), (0, 0), (0, 1)]  # (i, j)


def _frame() -> pd.DataFrame:
    return pd.DataFrame(
        {"x": X, "y": Y, "wt": np.arange(10.0), "z": np.arange(10.0), "g": ["a"] * 5 + ["b"] * 5}
    )


def _plot(df, layer, **aes_extra) -> p9.ggplot:
    return p9.ggplot(df, p9.aes("x", "y", **aes_extra)) + layer


def _first_rows(out: pd.DataFrame) -> pd.DataFrame:
    return out.drop_duplicates("group")


def _strict_warnings() -> None:
    # Fail on any warning except the deprecations pyproject already ignores:
    # pandas 3 raises Pandas4Warning from inside plotnine, not from ggann.
    warnings.simplefilter("error")
    warnings.simplefilter("ignore", DeprecationWarning)
    warnings.simplefilter("ignore", FutureWarning)


def _render(plot) -> None:
    with warnings.catch_warnings():
        _strict_warnings()
        plot.save(BytesIO(), format="png", verbose=False)


def test_exported_from_top_level_and_shapes():
    for name in ("geom_hex", "stat_binhex", "gghex"):
        assert getattr(ag, name) is getattr(shapes, name)
        assert name in ag.__all__
    assert not hasattr(ag, "bin_hexagons")


def test_bin_hexagons_table():
    out = bin_hexagons(X, Y, width=1, height=1, origin=(0, 0))
    assert out[["j", "i"]].values.tolist() == [[0, 0], [0, 1], [1, 0], [1, 1], [2, 2]]
    assert out["x"].tolist() == pytest.approx([0, 1, 0.5, 1.5, 2])
    assert out["y"].tolist() == pytest.approx([0, 0, 0.75, 0.75, 1.5])
    assert out["count"].tolist() == [3, 2, 3, 1, 1]
    assert out["count"].sum() == 10
    weighted = bin_hexagons(X, Y, width=1, height=1, origin=(0, 0), weight=np.arange(10))
    assert weighted["count"].tolist() == [9, 5, 18, 6, 7]
    z = bin_hexagons(X, Y, width=1, height=1, origin=(0, 0), z=np.arange(10), fun="mean")
    assert z["value"].tolist() == pytest.approx([3, 2.5, 6, 6, 7])


@pytest.mark.parametrize(
    ("fun", "expected"),
    [("sum", [9, 5, 18, 6, 7]), ("min", [0, 2, 4, 6, 7]), ("max", [8, 3, 9, 6, 7]),
     ("median", [1, 2.5, 5, 6, 7]), (lambda s: s.count(), [3, 2, 3, 1, 1])],
)  # fmt: skip
def test_bin_hexagons_summary_functions(fun, expected):
    out = bin_hexagons(X, Y, width=1, height=1, origin=(0, 0), z=np.arange(10), fun=fun)
    assert out["value"].tolist() == pytest.approx(expected)


def test_bin_hexagons_matches_brute_force_nearest_centre():
    rng = np.random.default_rng(1)
    x, y = rng.uniform(-3, 3, 2000), rng.uniform(-2, 2, 2000)
    w, h = 0.7, 0.9
    out = bin_hexagons(x, y, width=w, height=h, origin=(-3, -2))
    # Nearest centre in the metric where a regular hexagon is a Voronoi cell.
    ii, jj = np.meshgrid(np.arange(-6, 12), np.arange(-6, 12))
    cx = -3 + w * (ii + 0.5 * (jj % 2))
    cy = -2 + 0.75 * h * jj
    d = ((x[:, None] - cx.ravel()) / w) ** 2 + 0.75 * ((y[:, None] - cy.ravel()) / (0.75 * h)) ** 2
    nearest = d.argmin(axis=1)
    expected = pd.Series(list(zip(jj.ravel()[nearest], ii.ravel()[nearest]))).value_counts()
    got = {(j, i): c for j, i, c in zip(out["j"], out["i"], out["count"])}
    assert got == {k: v for k, v in expected.items()}


def test_tie_goes_to_the_lower_row():
    out = bin_hexagons([0.25], [0.375], width=1, height=1, origin=(0, 0))
    assert out[["j", "i"]].values.tolist() == [[0, 0]]


def test_hexagon_vertices_clockwise_from_top():
    vx, vy = hexagon_vertices(np.array([0.5]), np.array([0.75]), 1, 1)
    verts = list(zip(vx[0], vy[0]))
    assert verts == pytest.approx(
        [(0.5, 1.25), (1.0, 1.0), (1.0, 0.5), (0.5, 0.25), (0.0, 0.5), (0.0, 1.0)]
    )


def test_bin_hexagons_errors_and_empty():
    with pytest.raises(ValueError, match="finite"):
        bin_hexagons([np.nan], [0], width=1, height=1, origin=(0, 0))
    with pytest.raises(ValueError, match="positive"):
        bin_hexagons([0], [0], width=0, height=1, origin=(0, 0))
    with pytest.raises(ValueError, match="z needs fun"):
        bin_hexagons([0], [0], width=1, height=1, origin=(0, 0), z=[1])
    with pytest.raises(ValueError, match="fun must be"):
        bin_hexagons([0], [0], width=1, height=1, origin=(0, 0), z=[1], fun="mode")
    with pytest.raises(ValueError, match="z must be finite"):
        bin_hexagons([0, 1], [0, 1], width=1, height=1, origin=(0, 0), z=[1, np.nan], fun="mean")
    assert len(bin_hexagons([], [], width=1, height=1, origin=(0, 0))) == 0


def test_layer_table():
    out = _plot(_frame(), ag.geom_hex(binwidth=(1, 1))).layer_data(0)
    assert len(out) == 30
    assert out["group"].nunique() == 5
    first = _first_rows(out)
    assert first["count"].tolist() == [3, 2, 3, 1, 1]
    assert first["density"].tolist() == pytest.approx([0.3, 0.2, 0.3, 0.1, 0.1])
    assert first["ncount"].tolist() == pytest.approx([1, 2 / 3, 1, 1 / 3, 1 / 3])
    assert first[["i", "j"]].values.tolist() == [[0, 0], [1, 0], [0, 1], [1, 1], [2, 2]]
    assert out.loc[out["group"] == 3, "x"].tolist() == pytest.approx([0.5, 1, 1, 0.5, 0, 0])
    assert out.loc[out["group"] == 3, "y"].tolist() == pytest.approx([1.25, 1, 0.5, 0.25, 0.5, 1])
    assert "value" not in out.columns


def test_weight_and_z_aesthetics():
    weighted = _plot(_frame(), ag.geom_hex(binwidth=(1, 1)), weight="wt").layer_data(0)
    assert _first_rows(weighted)["count"].tolist() == [9, 5, 18, 6, 7]
    layer = ag.geom_hex(p9.aes(fill=p9.after_stat("value")), binwidth=(1, 1), fun="mean")
    out = _plot(_frame(), layer, z="z").layer_data(0)
    assert _first_rows(out)["value"].tolist() == pytest.approx([3, 2.5, 6, 6, 7])


def test_facets_share_one_lattice():
    plot = _plot(_frame(), ag.geom_hex(binwidth=(1, 1))) + p9.facet_wrap("g")
    out = plot.layer_data(0)
    panels = [_first_rows(g) for _, g in out.groupby("PANEL", observed=True)]
    assert [p["count"].tolist() for p in panels] == [[2, 2, 1], [1, 2, 1, 1]]
    assert all(g["x"].min() == pytest.approx(-0.5) for _, g in out.groupby("PANEL", observed=True))
    for p in panels:
        assert p["density"].sum() == pytest.approx(1)
    _render(plot)


def test_free_facet_scales_give_one_lattice_per_panel():
    plot = _plot(_frame(), ag.geom_hex(binwidth=(1, 1))) + p9.facet_wrap("g", scales="free")
    mins = [g["x"].min() for _, g in plot.layer_data(0).groupby("PANEL", observed=True)]
    assert mins[0] != pytest.approx(mins[1])


def test_missing_xy_dropped_with_warning():
    df = _frame()
    df.loc[0, "x"] = np.nan
    df.loc[1, "y"] = np.inf
    with pytest.warns(PlotnineWarning, match="Removed 2 rows"):
        out = _plot(df, ag.geom_hex(binwidth=(1, 1))).layer_data(0)
    assert _first_rows(out)["count"].sum() == 8
    with warnings.catch_warnings():
        _strict_warnings()
        _plot(df, ag.geom_hex(binwidth=(1, 1), na_rm=True)).layer_data(0)
    with pytest.raises(ValueError, match="finite"):
        bin_hexagons(df["x"], df["y"], width=1, height=1, origin=(0, 0))


def test_missing_z_raises():
    df = _frame()
    df.loc[3, "z"] = np.nan
    plot = _plot(df, ag.geom_hex(binwidth=(1, 1), fun="mean"), z="z")
    with pytest.raises(ValueError, match="z must be finite"):
        plot.layer_data(0)


def test_scale_limits_drop_points_and_move_origin():
    plot = _plot(_frame(), ag.geom_hex(binwidth=(1, 1))) + p9.scale_x_continuous(limits=(0.3, 3))
    with pytest.warns(PlotnineWarning):
        out = plot.layer_data(0)
    assert _first_rows(out)["count"].sum() == 7


def test_discrete_scale_raises_type_error():
    plot = p9.ggplot(_frame(), p9.aes("g", "y")) + ag.geom_hex(binwidth=(1, 1))
    with pytest.raises(TypeError, match="continuous"):
        plot.layer_data(0)


def test_degenerate_inputs():
    empty = _plot(_frame().iloc[:0], ag.geom_hex(bins=5))
    assert len(empty.layer_data(0)) == 0
    _render(empty)
    one = _plot(_frame().iloc[:1], ag.geom_hex(bins=5)).layer_data(0)
    assert one["group"].nunique() == 1
    same = pd.DataFrame({"x": [1.0, 1.0], "y": [2.0, 2.0]})
    out = _plot(same, ag.geom_hex(bins=5)).layer_data(0)
    assert out["group"].nunique() == 1
    assert out["count"].iloc[0] == 2


def test_scalar_bins_are_regular_hexagons():
    df = pd.DataFrame({"x": [0.0, 10.0], "y": [0.0, 3.0]})
    out = _plot(df, ag.geom_hex(bins=10)).layer_data(0)
    width = out["x"].max() - out["x"].min()
    height = out["y"].max() - out["y"].min()
    one = out[out["group"] == 1]
    assert one["x"].max() - one["x"].min() == pytest.approx(1.0)
    assert one["y"].max() - one["y"].min() == pytest.approx(2 / np.sqrt(3))
    assert width > 0 and height > 0
    pair = _plot(df, ag.geom_hex(bins=(10, 3))).layer_data(0)
    one = pair[pair["group"] == 1]
    assert one["y"].max() - one["y"].min() == pytest.approx(1.0)
    wins = _plot(df, ag.geom_hex(bins=3, binwidth=2)).layer_data(0)
    one = wins[wins["group"] == 1]
    assert one["x"].max() - one["x"].min() == pytest.approx(2.0)


def test_seam_rule():
    out = _plot(_frame(), ag.geom_hex(binwidth=(1, 1))).layer_data(0)
    assert (out["color"] == out["fill"]).all()
    assert (out["size"] == 0.2).all()
    red = _plot(_frame(), ag.geom_hex(binwidth=(1, 1), fill="red")).layer_data(0)
    assert (red["color"] == "red").all()
    half = _plot(_frame(), ag.geom_hex(binwidth=(1, 1), alpha=0.5)).layer_data(0)
    assert half["color"].isna().all()
    kept = _plot(_frame(), ag.geom_hex(binwidth=(1, 1), color="black")).layer_data(0)
    assert (kept["color"] == "black").all()
    kept = _plot(_frame(), ag.geom_hex(p9.aes(colour="x"), binwidth=(1, 1))).layer_data(0)
    assert not (kept["color"] == kept["fill"]).all()


def test_missing_and_bad_sizes_raise_at_build():
    with pytest.raises(ValueError, match="no default"):
        _plot(_frame(), ag.geom_hex()).layer_data(0)
    for kwargs in ({"bins": 0}, {"binwidth": (0, 1)}, {"bins": 2.5}, {"bins": (3, 3, 3)}):
        with pytest.raises(ValueError):
            _plot(_frame(), ag.geom_hex(**kwargs)).layer_data(0)
    with pytest.raises(ValueError, match="fun"):
        _plot(_frame(), ag.geom_hex(bins=3), z="z").layer_data(0)
    with pytest.raises(ValueError, match="fun must be"):
        _plot(_frame(), ag.geom_hex(bins=3, fun="mode"), z="z").layer_data(0)
    with pytest.raises(PlotnineError):
        _plot(_frame(), ag.geom_hex(bin=3))


def test_renders_under_coord_facet_and_log_scales():
    base = _plot(_frame().assign(x=lambda d: d["x"] + 1), ag.geom_hex(bins=4))
    _render(base + p9.coord_equal())
    _render(base + p9.coord_flip())
    _render(base + p9.scale_x_log10())
    _render(base + p9.facet_wrap("g"))


def test_gghex_returns_a_plot_with_counts():
    plot = ag.gghex(_frame(), "x", "y", binwidth=(1, 1))
    assert isinstance(plot, p9.ggplot)
    assert isinstance(plot.coordinates, p9.coord_equal)
    assert _first_rows(plot.layer_data(0))["count"].tolist() == [3, 2, 3, 1, 1]
    _render(plot)
    assert not isinstance(
        ag.gghex(_frame(), "x", "y", bins=4, equal_aspect=False).coordinates, p9.coord_equal
    )


def test_gghex_stat_weight_z_and_facet():
    plot = ag.gghex(_frame(), "x", "y", binwidth=(1, 1), weight="wt", stat="ncount")
    assert _first_rows(plot.layer_data(0))["fill"].nunique() == 5
    plot = ag.gghex(_frame(), "x", "y", z="z", fun="mean", binwidth=(1, 1), stat="value", facet="g")
    out = plot.layer_data(0)
    assert out["PANEL"].nunique() == 2
    assert "value" in out.columns
    _render(plot)
    call = ag.gghex(_frame(), "x", "y", z="z", fun=np.max, bins=3, stat="value")
    assert _first_rows(call.layer_data(0))["value"].max() == 9


def test_gghex_errors():
    df = _frame()
    with pytest.raises(KeyError):
        ag.gghex(df, "nope", "y", bins=3)
    with pytest.raises(KeyError):
        ag.gghex(df, "x", "y", bins=3, facet="nope")
    with pytest.raises(ValueError, match="stat must be"):
        ag.gghex(df, "x", "y", bins=3, stat="median")
    with pytest.raises(ValueError, match="needs z"):
        ag.gghex(df, "x", "y", bins=3, stat="value")
    with pytest.raises(ValueError, match="needs fun"):
        ag.gghex(df, "x", "y", bins=3, z="z")
    with pytest.raises(ValueError, match="fun must be"):
        ag.gghex(df, "x", "y", bins=3, z="z", fun="mode")
    with pytest.raises(ValueError, match="no default"):
        ag.gghex(df, "x", "y")
    with pytest.raises(ValueError):
        ag.gghex(df, "x", "y", bins=0)
    with pytest.raises(ValueError):
        ag.gghex(df, "x", "y", binwidth=(0, 1))


# plot_embedding ---------------------------------------------------------------


def test_embedding_hex_counts_cells(adata):
    plot = ag.plot_embedding(adata, "umap", hex=True, hex_bins=12)
    assert isinstance(plot, p9.ggplot)
    assert len(plot.data) == adata.n_obs
    out = plot.layer_data(0)
    assert _first_rows(out)["count"].sum() == adata.n_obs
    assert plot.labels.fill == "cells"
    _render(plot)


def test_embedding_hex_split_by_counts_per_panel(adata, group_key):
    plot = ag.plot_embedding(adata, "umap", hex=True, hex_bins=10, split_by=group_key)
    out = plot.layer_data(0)
    per_panel = [_first_rows(g)["count"].sum() for _, g in out.groupby("PANEL", observed=True)]
    sizes = adata.obs[group_key].value_counts()
    assert sorted(per_panel) == sorted(sizes[sizes > 0].tolist())
    # Fixed scales: every hexagon centre sits on one lattice anchored at the x minimum.
    one = out.groupby(["PANEL", "group"], observed=True).get_group((out["PANEL"].iloc[0], 1))
    width = one["x"].max() - one["x"].min()
    centres = out.groupby(["PANEL", "group"], observed=True)["x"].first()
    steps = (centres - adata.obsm["X_umap"][:, 0].min()) / width * 2
    assert np.allclose(steps, steps.round(), atol=1e-6)
    _render(plot)


def test_embedding_hex_numeric_color(adata, markers):
    plot = ag.plot_embedding(adata, "umap", markers[0], hex=True, hex_bins=10, hex_fun="mean")
    out = plot.layer_data(0)
    assert "value" in out.columns
    assert plot.labels.fill == markers[0]
    _render(plot)
    sums = ag.plot_embedding(adata, "umap", "n_counts", hex=True, hex_bins=10, hex_fun="sum")
    assert _first_rows(sums.layer_data(0))["value"].sum() == pytest.approx(
        adata.obs["n_counts"].sum()
    )
    _render(
        ag.plot_embedding(
            adata, "umap", "n_counts", hex=True, hex_bins=10, hex_fun="max", color_map="magma"
        )
    )


def test_embedding_hex_color_list_forwards_options(adata, markers):
    composed = ag.plot_embedding(adata, "umap", markers[:2], hex=True, hex_bins=8, hex_fun="median")
    _render(composed)
    layers = [p.layer_data(0) for p in composed]
    assert len(layers) == 2
    assert all("value" in out.columns for out in layers)


def test_embedding_hex_rejections(adata, group_key, markers):
    kw = {"hex": True, "hex_bins": 8}
    with pytest.raises(ValueError, match="matplotlib"):
        ag.plot_embedding(adata, "umap", backend="matplotlib", **kw)
    with pytest.raises(ValueError, match="hex_bins"):
        ag.plot_embedding(adata, "umap", hex=True)
    with pytest.raises(ValueError, match="hex_bins"):
        ag.plot_embedding(adata, "umap", hex=True, hex_bins=0)
    with pytest.raises(ValueError, match="categorical"):
        ag.plot_embedding(adata, "umap", group_key, hex_fun="mean", **kw)
    with pytest.raises(ValueError, match="hex_fun"):
        ag.plot_embedding(adata, "umap", markers[0], **kw)
    with pytest.raises(ValueError, match="fun must be"):
        ag.plot_embedding(adata, "umap", markers[0], hex_fun="mode", **kw)
    with pytest.raises(ValueError, match="downsample"):
        ag.plot_embedding(adata, "umap", downsample=50, **kw)
    with pytest.raises(ValueError, match="vmin"):
        ag.plot_embedding(adata, "umap", markers[0], hex_fun="mean", vmax=2, **kw)
    with pytest.raises(ValueError, match="method"):
        ag.plot_embedding(adata, "umap", method="binned", **kw)
    with pytest.raises(ValueError, match="pointdensity"):
        ag.plot_embedding(adata, "umap", pointdensity=False, **kw)
    with pytest.raises(ValueError, match="hex=True"):
        ag.plot_embedding(adata, "umap", hex_bins=8)
    with pytest.raises(ValueError, match="hex=True"):
        ag.plot_embedding(adata, "umap", hex_fun="mean")


def test_embedding_hex_missing_numeric_value_raises(adata):
    sub = adata.copy()
    sub.obs["score"] = np.linspace(0, 1, sub.n_obs)
    sub.obs.loc[sub.obs.index[0], "score"] = np.nan
    with pytest.raises(ValueError, match="finite"):
        ag.plot_embedding(sub, "umap", "score", hex=True, hex_bins=8, hex_fun="mean")
