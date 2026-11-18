"""stat_pie, stat_radar, ggpie, ggring, ggradar and plot_proportions kinds."""

from __future__ import annotations

import warnings
from io import BytesIO

import numpy as np
import pandas as pd
import plotnine as p9
import pytest

import ggann as ag
import ggann.shapes as shapes


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


def _pie(weights, **kwargs):
    df = pd.DataFrame({"g": list("abcdefgh")[: len(weights)], "w": weights})
    return p9.ggplot(df, p9.aes(fill="g", weight="w")) + ag.stat_pie(**kwargs)


def _radar_frame() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "axis": list("abcde") * 2,
            "v": [4, 2, 3, 1, 2, 1, 3, 2, 4, 3],
            "s": ["g1"] * 5 + ["g2"] * 5,
        }
    )


def test_exported_from_top_level_and_shapes():
    names = [
        "geom_pie",
        "stat_pie",
        "geom_radar",
        "stat_radar",
        "coord_radar",
        "ggpie",
        "ggring",
        "ggradar",
    ]
    for name in names:
        assert getattr(ag, name) is getattr(shapes, name)
        assert name in ag.__all__


# stat_pie -----------------------------------------------------------------------


def test_pie_table_for_weights():
    out = _pie([4, 3, 2, 1]).layer_data(0)
    assert out["count"].tolist() == pytest.approx([4, 3, 2, 1])
    assert out["prop"].tolist() == pytest.approx([0.4, 0.3, 0.2, 0.1])
    assert out["ymin"].tolist() == pytest.approx([0, 0.4, 0.7, 0.9])
    assert out["ymax"].tolist() == [pytest.approx(v) for v in (0.4, 0.7, 0.9, 1.0)]
    assert out["y"].tolist() == pytest.approx([0.2, 0.55, 0.8, 0.95])
    assert out["xmin"].tolist() == pytest.approx([0.5] * 4)
    assert out["xmax"].tolist() == pytest.approx([1.5] * 4)
    assert out["x"].tolist() == pytest.approx([1.0] * 4)
    assert out["ymax"].iloc[-1] == 1.0


def test_pie_counts_raw_rows():
    df = pd.DataFrame({"g": list("aaabbc")})
    out = (p9.ggplot(df, p9.aes(fill="g")) + ag.stat_pie()).layer_data(0)
    assert out["count"].tolist() == pytest.approx([3, 2, 1])
    assert out["prop"].tolist() == pytest.approx([1 / 2, 1 / 3, 1 / 6])
    assert out["ymax"].tolist() == pytest.approx([0.5, 5 / 6, 1.0])


def test_pie_closes_in_each_facet():
    df = pd.DataFrame({"g": list("abc") * 2, "s": list("xxxyyy"), "w": [1, 1, 2, 5, 3, 2]})
    plot = p9.ggplot(df, p9.aes(fill="g", weight="w")) + ag.stat_pie() + p9.facet_wrap("~s")
    out = plot.layer_data(0)
    first, second = (out[out["PANEL"] == k] for k in out["PANEL"].unique())
    assert first["ymax"].tolist() == pytest.approx([0.25, 0.5, 1.0])
    assert second["ymax"].tolist() == pytest.approx([0.5, 0.8, 1.0])


def test_two_rings_with_width():
    df = pd.DataFrame({"x": [1, 1, 2, 2], "g": list("abab"), "w": [1, 1, 3, 1]})
    out = (p9.ggplot(df, p9.aes("x", fill="g", weight="w")) + ag.stat_pie(width=0.9)).layer_data(0)
    assert out["xmin"].tolist() == pytest.approx([0.55, 0.55, 1.55, 1.55])
    assert out["xmax"].tolist() == pytest.approx([1.45, 1.45, 2.45, 2.45])
    assert out["ymax"].tolist() == pytest.approx([0.5, 1.0, 0.75, 1.0])


def test_zero_slices_stay_unless_dropped():
    out = _pie([4, 0, 2, 0]).layer_data(0)
    assert out["ymax"].tolist() == pytest.approx([2 / 3, 2 / 3, 1.0, 1.0])
    kept = _pie([4, 0, 2, 0], drop_zero=True).layer_data(0)
    assert len(kept) == 2


def test_single_category_is_a_full_circle():
    df = pd.DataFrame({"g": ["a"]})
    plot = (
        p9.ggplot(df, p9.aes(fill="g")) + ag.geom_pie() + p9.coord_radial(theta="y", expand=False)
    )
    out = plot.layer_data(0)
    assert out[["ymin", "ymax"]].iloc[0].tolist() == [0, 1]
    _render(plot)


def test_negative_weight_raises():
    with pytest.raises(ValueError, match="negative"):
        _pie([4, -1, 2]).layer_data(0)


def test_nan_weight_is_removed_with_a_warning():
    with pytest.warns(UserWarning, match="non-finite"):
        out = _pie([4, np.nan, 2, 1]).layer_data(0)
    assert out["count"].tolist() == pytest.approx([4, 2, 1])
    with warnings.catch_warnings():
        _strict_warnings()
        quiet = _pie([4, np.nan, 2, 1], na_rm=True).layer_data(0)
    assert quiet["count"].tolist() == pytest.approx([4, 2, 1])


def test_ring_with_zero_total_is_dropped_with_a_warning():
    with pytest.warns(UserWarning, match="total weight 0"):
        out = _pie([0, 0]).layer_data(0)
    assert len(out) == 0


def test_label_position_moves_the_label_radius():
    out = _pie([1, 1], label_position=0.25).layer_data(0)
    assert out["x"].tolist() == pytest.approx([0.75, 0.75])


# stat_radar ---------------------------------------------------------------------


def _radar_layer(df, **kwargs):
    plot = p9.ggplot(df, p9.aes("axis", "v", group="s")) + ag.stat_radar(**kwargs)
    return plot.layer_data(0)


def test_radar_sorts_and_closes():
    df = _radar_frame().iloc[:5].copy()
    out = _radar_layer(df)
    assert len(out) == 6
    assert out["x"].tolist() == [1, 2, 3, 4, 5, 1]
    assert out["y"].tolist() == [4, 2, 3, 1, 2, 4]
    shuffled = _radar_layer(df.sample(frac=1, random_state=1))
    pd.testing.assert_frame_equal(shuffled.reset_index(drop=True), out.reset_index(drop=True))
    assert len(_radar_layer(df, closed=False)) == 5


# ggpie and ggring ---------------------------------------------------------------


def test_ggpie_and_ggring_render():
    df = pd.DataFrame({"k": list("abcd"), "n": [4, 3, 2, 1], "s": list("xyxy")})
    _render(ag.ggpie(df, "k", "n"))
    _render(ag.ggring(df, "k", "n", inner_radius=0.6, label="count"))
    _render(ag.ggpie(df, "k", split_by="s"))


def test_ggpie_percent_label_skips_a_zero_slice():
    df = pd.DataFrame({"k": list("abcd"), "n": [4, 0, 2, 0]})
    plot = ag.ggpie(df, "k", "n", label="percent")
    assert plot.layer_data(1)["label"].tolist() == ["67%", "33%"]
    _render(plot)


def test_ggpie_count_label_has_no_decimals():
    df = pd.DataFrame({"k": list("abc")})
    labels = ag.ggpie(df, "k", label="count").layer_data(1)["label"]
    assert labels.tolist() == ["1", "1", "1"]


def test_ggpie_wedge_angles_and_hole():
    df = pd.DataFrame({"k": list("abcd"), "n": [4, 3, 2, 1]})
    ax = ag.ggpie(df, "k", "n", inner_radius=0.6).draw().axes[0]
    vertices = ax.collections[0].get_paths()[0].vertices
    assert np.ptp(vertices[:, 0]) / (2 * np.pi) == pytest.approx(0.4, abs=1e-3)
    assert ax.get_rorigin() == pytest.approx(-1.0)


def test_ggpie_order_and_palette():
    df = pd.DataFrame({"k": list("abc"), "n": [1, 1, 1]})
    plot = ag.ggpie(
        df, "k", "n", order=["c", "a", "b"], palette={"a": "red", "b": "blue", "c": "green"}
    )
    assert plot.layer_data(0)["fill"].tolist() == ["green", "red", "blue"]
    listed = ag.ggpie(df, "k", "n", palette=["#111111", "#222222", "#333333"])
    assert listed.layer_data(0)["fill"].tolist() == ["#111111", "#222222", "#333333"]


def test_ggpie_drops_missing_category_with_a_warning():
    df = pd.DataFrame({"k": ["a", None, "b"], "n": [1, 1, 1]})
    with pytest.warns(UserWarning, match="dropped 1 of 3 rows"):
        plot = ag.ggpie(df, "k", "n")
    assert len(plot.layer_data(0)) == 2


@pytest.mark.parametrize(
    "call, match",
    [
        (lambda df: ag.ggring(df, "k", inner_radius=0), "inner_radius"),
        (lambda df: ag.ggpie(df, "k", inner_radius=1), "inner_radius"),
        (lambda df: ag.ggpie(df, "k", label="bad"), "label"),
        (lambda df: ag.ggpie(df, "k", palette={"a": "red"}), "no colour"),
        (lambda df: ag.ggpie(df, "k", start=np.nan), "start"),
    ],
)
def test_ggpie_errors(call, match):
    with pytest.raises(ValueError, match=match):
        call(pd.DataFrame({"k": list("ab")}))


def test_ggpie_accepts_the_scanpy_spelling():
    df = pd.DataFrame({"k": list("ab")})
    assert len(ag.ggpie(df, groupby="k").layer_data(0)) == 2


# ggradar ------------------------------------------------------------------------


def test_ggradar_geometry():
    ax = ag.ggradar(_radar_frame(), "axis", "v", "s").draw().axes[0]
    assert ax.get_rorigin() == 0
    assert ax.get_xticks().tolist() == pytest.approx([2 * np.pi * i / 5 for i in range(5)])
    labels = [t.get_text() for t in ax.get_yticklabels()]
    assert labels == [f"{t:g}" for t in ax.get_yticks()]
    assert [len(p.vertices) for p in ax.collections[0].get_paths()] == [7, 7]


def test_ggradar_single_series_and_palette():
    df = _radar_frame().iloc[:5]
    _render(ag.ggradar(df, "axis", "v"))
    plot = ag.ggradar(_radar_frame(), "axis", "v", "s", palette={"g1": "red", "g2": "blue"})
    assert plot.layer_data(0)["fill"].unique().tolist() == ["red", "blue"]
    _render(plot)


def test_ggradar_rlim_and_order():
    plot = ag.ggradar(
        _radar_frame(), "axis", "v", "s", rlim=(0, 6), axis_order=list("edcba"), points=False
    )
    assert plot.layer_data(0)["x"].iloc[:5].tolist() == [1, 2, 3, 4, 5]
    assert plot.layer_data(0)["y"].iloc[:5].tolist() == [2, 1, 3, 2, 4]
    _render(plot)


@pytest.mark.parametrize(
    "edit, match",
    [
        (lambda d: d.iloc[:-1], "pair needs a row"),
        (lambda d: pd.concat([d, d.iloc[:1]]), "more than one row"),
        (lambda d: d[d["axis"].isin(["a", "b"])], "at least 3"),
        (lambda d: d.assign(v=d["v"].mask(d.index == 0)), "missing or non-finite"),
        (lambda d: d.assign(v=d["v"].where(d.index != 0, -1.0)), "outside rlim"),
        (lambda d: d.assign(axis=d["axis"].mask(d.index == 0)), "missing labels"),
    ],
)
def test_ggradar_errors(edit, match):
    with pytest.raises(ValueError, match=match):
        ag.ggradar(edit(_radar_frame().reset_index(drop=True)), "axis", "v", "s")


def test_ggradar_rejects_a_value_above_rlim_and_a_bad_rlim():
    with pytest.raises(ValueError, match="outside rlim"):
        ag.ggradar(_radar_frame(), "axis", "v", "s", rlim=(0, 3))
    with pytest.raises(ValueError, match="increasing"):
        ag.ggradar(_radar_frame(), "axis", "v", "s", rlim=(2, 2))


# plot_proportions ---------------------------------------------------------------


@pytest.mark.parametrize("kind", ["pie", "ring"])
@pytest.mark.parametrize("split", [None, "depth"])
def test_proportions_pie_and_ring_build(adata, group_key, kind, split):
    base = ag.plot_proportions(adata, group_key, split, kind="bar")
    plot = ag.plot_proportions(adata, group_key, split, kind=kind)
    assert list(plot.data.columns) == list(base.data.columns)
    out = plot.layer_data(0)
    assert out.groupby("PANEL", observed=True)["ymax"].max().tolist() == pytest.approx(
        [1.0] * (2 if split else 1)
    )
    _render(plot)


@pytest.mark.parametrize("split", [None, "depth"])
def test_proportions_radar_builds(adata, group_key, split):
    kwargs = {"radar_norm": "within_split"} if split else {}
    plot = ag.plot_proportions(adata, group_key, split, kind="radar", **kwargs)
    base = ag.plot_proportions(adata, group_key, split)
    assert list(plot.data.columns) == list(base.data.columns)
    _render(plot)


def test_proportions_radar_norm_changes_the_numbers(adata, group_key):
    within = ag.plot_proportions(
        adata, group_key, "depth", kind="radar", radar_norm="within_split"
    ).data
    across = ag.plot_proportions(
        adata, group_key, "depth", kind="radar", radar_norm="across_splits"
    ).data
    assert within.groupby("depth", observed=True)["value"].sum().tolist() == pytest.approx([1, 1])
    assert across.groupby(group_key, observed=True)["value"].sum().tolist() == pytest.approx(
        [1] * across[group_key].nunique()
    )
    assert not np.allclose(within["value"], across["value"])


def test_proportions_radar_needs_a_normalisation_choice(adata, group_key):
    with pytest.raises(ValueError, match="radar_norm"):
        ag.plot_proportions(adata, group_key, "depth", kind="radar")
    with pytest.raises(ValueError, match="radar_norm"):
        ag.plot_proportions(adata, group_key, "depth", kind="radar", radar_norm="bad")
    counts = ag.plot_proportions(adata, group_key, "depth", kind="radar", normalize=False)
    assert counts.data["value"].equals(counts.data["n"])


def test_proportions_radar_needs_three_groups(adata):
    two = adata[adata.obs["bulk_labels"].isin(["CD56+ NK", "CD14+ Monocyte"])].copy()
    with pytest.raises(ValueError, match="at least 3"):
        ag.plot_proportions(two, "bulk_labels", kind="radar")


def test_proportions_unknown_kind_lists_the_options(adata, group_key):
    with pytest.raises(ValueError, match="pie.*ring.*radar"):
        ag.plot_proportions(adata, group_key, kind="waffle")
