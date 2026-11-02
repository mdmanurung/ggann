"""Regressions for the 0.1.0 readiness rows owned by the ``dist`` area.

Row IDs refer to plans/evaluation/remaining_todos.md.
"""

from __future__ import annotations

import hashlib
import io
import warnings

import numpy as np
import pandas as pd
import pytest
from plotnine import geom_col, labs, scale_y_continuous, theme
from plotnine.exceptions import PlotnineWarning

import ggann as ag

GROUP = "bulk_labels"
GENES = ["CD3D", "NKG7"]


def _png(plot) -> str:
    import matplotlib.pyplot as plt

    buf = io.BytesIO()
    plot.save(buf, format="png", dpi=40, verbose=False)
    plt.close("all")
    return hashlib.md5(buf.getvalue()).hexdigest()


def _scatter_counts(plot) -> list[int]:
    """Number of marks in each point collection per axes (violin polygons excluded)."""
    import matplotlib.pyplot as plt

    fig = plot.draw()
    out = [
        len(c.get_offsets()) for ax in fig.axes for c in ax.collections if len(c.get_offsets()) > 1
    ]
    plt.close(fig)
    return out


def _key(adata) -> tuple:
    x = adata.X if isinstance(adata.X, np.ndarray) else adata.X.toarray()
    obs = pd.util.hash_pandas_object(adata.obs, index=True).to_numpy()
    return hashlib.sha256(x.tobytes() + obs.tobytes()).hexdigest(), list(adata.obs.columns)


# --------------------------------------------------------------------------- #
# B6c: ridge n_grid / scale validation
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("bad", [0, 1, -3, 2.5, True, "a", None])
def test_B6c_ridge_n_grid_rejected(adata, bad):
    with pytest.raises(ValueError, match="n_grid"):
        ag.plot_ridge(adata, GENES[:1], GROUP, use_raw=True, n_grid=bad)


@pytest.mark.parametrize("bad", [0, -1.0, float("nan"), float("inf"), True, "x", None])
def test_B6c_ridge_scale_rejected(adata, bad):
    with pytest.raises(ValueError, match="scale"):
        ag.plot_ridge(adata, GENES[:1], GROUP, use_raw=True, scale=bad)


def test_B6c_ridge_valid_still_builds(adata):
    ag.plot_ridge(adata, GENES[:1], GROUP, use_raw=True, n_grid=2, scale=0.5).draw()


# --------------------------------------------------------------------------- #
# A4: numeric obs keys, typed selectors, log1p, order, rotation, stripplot
# --------------------------------------------------------------------------- #
OBS_CALLERS = {
    "violin": lambda a, f: ag.plot_violin(a, f, GROUP),
    "box": lambda a, f: ag.plot_box(a, f, GROUP),
    "sina": lambda a, f: ag.plot_sina(a, f, GROUP),
    "ridge": lambda a, f: ag.plot_ridge(a, f, GROUP, n_grid=16),
    "stacked": lambda a, f: ag.plot_stacked_violin(a, f, GROUP),
    "bar": lambda a, f: ag.plot_expression_bar(a, f, GROUP),
}


@pytest.mark.parametrize("name", OBS_CALLERS)
@pytest.mark.parametrize("sel", ["bare", "typed", "prefixed"])
def test_A4_numeric_obs_key_is_plotted(adata, name, sel):
    spec = {"bare": "n_genes", "typed": ag.obs("n_genes"), "prefixed": "obs:n_genes"}[sel]
    plot = OBS_CALLERS[name](adata, [spec])
    plot.draw()
    if name in ("ridge", "bar"):  # plotted frames are summaries; the feature is named
        assert list(plot.data["feature"].astype(str).unique()) == ["n_genes"]
        return
    data = plot.data
    assert list(data["feature"].astype(str).unique()) == ["n_genes"]
    assert np.array_equal(data["value"].to_numpy(), adata.obs["n_genes"].to_numpy(dtype=float))
    assert (data["obs_name"].to_numpy() == adata.obs_names.to_numpy()).all()


def test_A4_gene_and_obs_mix_keeps_order_and_alignment(adata):
    plot = ag.plot_violin(adata, ["n_genes", "CD3D"], GROUP, use_raw=True)
    d = plot.data
    assert list(d["feature"].cat.categories) == ["n_genes", "CD3D"]
    cd3d = d[d["feature"] == "CD3D"].set_index("obs_name")["value"]
    expect = adata.raw[:, "CD3D"].X
    expect = expect.toarray().ravel() if hasattr(expect, "toarray") else np.ravel(expect)
    assert np.allclose(cd3d.reindex(adata.obs_names).to_numpy(), expect)
    ngenes = d[d["feature"] == "n_genes"].set_index("obs_name")["value"]
    assert np.allclose(ngenes.reindex(adata.obs_names).to_numpy(), adata.obs["n_genes"])


def test_A4_bare_name_in_both_stays_a_gene(adata):
    a = adata.copy()
    a.obs["CD3D"] = 1.0  # collides with a gene name
    d = ag.plot_violin(a, ["CD3D"], GROUP, use_raw=True).data
    assert d["value"].nunique() > 1  # the gene, not the constant obs column
    d = ag.plot_violin(a, [ag.obs("CD3D")], GROUP, use_raw=True).data
    assert set(d["value"]) == {1.0}


def test_A4_bad_selectors_fail_clearly(adata):
    with pytest.raises(TypeError, match="not numeric"):
        ag.plot_violin(adata, [ag.obs(GROUP)], GROUP)
    with pytest.raises(KeyError, match="no_such_col"):
        ag.plot_violin(adata, [ag.obs("no_such_col")], GROUP)
    with pytest.raises(ValueError, match="per-gene layer"):
        ag.plot_violin(adata, [ag.gene("CD3D", layer="x")], GROUP)
    with pytest.raises(ValueError, match="duplicate"):
        ag.plot_violin(adata, ["n_genes", ag.obs("n_genes")], GROUP)


def test_A4_log_is_log1p_with_zeros_kept(adata):
    raw = ag.plot_violin(adata, GENES, GROUP, use_raw=True).data
    logged = ag.plot_violin(adata, GENES, GROUP, use_raw=True, log=True)
    assert (raw["value"] == 0).any()
    assert np.allclose(logged.data["value"], np.log1p(raw["value"]))
    assert (logged.data["value"] == 0).sum() == (raw["value"] == 0).sum()
    assert logged.labels.y == "log1p(expression)"
    for fn in (ag.plot_box, ag.plot_sina):
        assert fn(adata, GENES, GROUP, use_raw=True, log=True).labels.y == "log1p(expression)"


def test_A4_log_rejects_scaled_negative_values(adata):
    a = adata.copy()
    a.obs["neg"] = np.linspace(-3, 3, a.n_obs)
    with pytest.raises(ValueError, match="log1p"):
        ag.plot_violin(a, [ag.obs("neg")], GROUP, log=True)


def test_A4_stacked_violin_scanpy_orientation_and_median_fill(adata):
    # default = Scanpy's default: groups on the facet rows, genes along x
    plot = ag.plot_stacked_violin(adata, GENES, GROUP, use_raw=True)
    assert type(plot.facet).__name__ == "facet_grid" and plot.facet.rows == [GROUP]
    d = plot.data
    med = d.groupby(["feature", GROUP], observed=True)["value"].median()
    got = d.groupby(["feature", GROUP], observed=True)["median"].first()
    assert np.allclose(med.to_numpy(), got.to_numpy())
    assert plot.mapping["x"] == "feature" and plot.mapping["fill"] == "median"
    fig = plot.draw()
    assert len(fig.axes) >= adata.obs[GROUP].nunique()  # one row per group


def test_stacked_violin_swap_axes_is_the_genes_as_rows_layout(adata):
    plot = ag.plot_stacked_violin(adata, GENES, GROUP, use_raw=True, swap_axes=True)
    assert plot.facet.rows == ["feature"]
    assert plot.mapping["x"] == GROUP and plot.mapping["fill"] == "median"
    same = ag.plot_stacked_violin(adata, GENES, GROUP, use_raw=True).data
    pd.testing.assert_frame_equal(plot.data, same)  # only the layout differs
    fig = plot.draw()
    assert len(fig.axes) >= len(GENES)  # one row per gene


def test_A4_order_subset_rotation_stripplot(adata):
    groups = list(adata.obs[GROUP].cat.categories[:3])[::-1]
    for fn in (ag.plot_violin, ag.plot_box, ag.plot_sina, ag.plot_stacked_violin):
        plot = fn(adata, GENES, GROUP, use_raw=True, order=groups, rotation=90)
        assert list(plot.data[GROUP].cat.categories) == groups
        plot.draw()
    with pytest.raises(ValueError, match="unknown"):
        ag.plot_violin(adata, GENES, GROUP, order=["nope"])
    with pytest.raises(ValueError, match="not both"):
        ag.plot_violin(adata, GENES, GROUP, order=groups, categories_order=groups)
    plain = ag.plot_stacked_violin(adata, GENES, GROUP, use_raw=True)
    strip = ag.plot_stacked_violin(adata, GENES, GROUP, use_raw=True, stripplot=True)
    assert len(strip.layers) == len(plain.layers) + 1
    assert sum(_scatter_counts(strip)) == adata.n_obs * len(GENES)


# --------------------------------------------------------------------------- #
# A3: a gene listed under two marker groups appears once per occurrence
# --------------------------------------------------------------------------- #
REPEAT = {"T": ["CD3D", "CD8A", "NKG7"], "NK": ["NKG7", "GNLY", "CD3D"]}


@pytest.mark.parametrize(
    ("grouped", "single"),
    [
        (ag.plot_dotplot_grouped, ag.plot_dotplot),
        (ag.plot_matrixplot_grouped, ag.plot_matrixplot),
    ],
)
def test_A3_repeated_gene_once_per_occurrence_values_from_one_computation(adata, grouped, single):
    plot = grouped(adata, REPEAT, GROUP, use_raw=True)
    d = plot.data
    per_panel = d.groupby("gene_group", observed=True)["gene"].apply(
        lambda s: list(dict.fromkeys(s))
    )
    assert per_panel["T"] == REPEAT["T"]  # per-panel order as listed
    assert per_panel["NK"] == REPEAT["NK"]
    once = single(adata, ["CD3D", "CD8A", "NKG7", "GNLY"], GROUP, use_raw=True).data
    col = "mean_expression"
    for g in ("CD3D", "NKG7"):
        ref = once[once["feature"] == g].set_index(GROUP)[col]
        for label in ("T", "NK"):
            part = d[(d["gene"] == g) & (d["gene_group"] == label)].set_index(GROUP)[col]
            assert np.allclose(part.reindex(ref.index).to_numpy(), ref.to_numpy())
    plot.draw()


def test_A3_tracksplot_dict_repeats_tracks(adata):
    plot = ag.plot_tracksplot(adata, REPEAT, GROUP, use_raw=True)
    tracks = list(plot.data["feature"].cat.categories)
    assert tracks == [f"{g} [{k}]" for k, gs in REPEAT.items() for g in gs]
    one = ag.plot_tracksplot(adata, ["NKG7"], GROUP, use_raw=True).data["value"].to_numpy()
    for track in ("NKG7 [T]", "NKG7 [NK]"):
        part = plot.data[plot.data["feature"] == track]["value"].to_numpy()
        assert np.array_equal(np.sort(part), np.sort(one))
    plot.draw()


def test_A3_grouped_rejects_non_dict(adata):
    with pytest.raises(TypeError, match="dict"):
        ag.plot_dotplot_grouped(adata, ["CD3D"], GROUP)


# --------------------------------------------------------------------------- #
# P1: tracksplot rectangles equal the old geom_col drawing
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("use_raw", [True, False])  # False: scaled, negative baselines
def test_P1_tracksplot_rects_match_geom_col(adata, use_raw):
    import matplotlib.pyplot as plt
    from plotnine import aes, facet_grid, ggplot

    plot = ag.plot_tracksplot(adata, GENES, GROUP, use_raw=use_raw)
    assert type(plot.layers[0].geom).__name__ == "geom_rect"
    data = plot.data
    old = (
        ggplot(data, aes("cell_rank", "value", fill=GROUP))
        + geom_col(width=1.0)
        + facet_grid("feature ~ .", scales="free_y")
        + ag.scale_fill_obs(adata, GROUP)
    )

    def rects(p):
        fig = p.draw()
        out = []
        for ax in fig.axes:
            for c in ax.collections:
                ext = np.array([pp.get_extents().extents for pp in c.get_paths()])
                fc = np.asarray(c.get_facecolor())
                fc = np.repeat(fc, len(ext), axis=0) if len(fc) == 1 else fc
                t = np.hstack([ext, fc])
                out.append(t[np.lexsort(t.T[::-1])])
        plt.close(fig)
        return out

    new_r, old_r = rects(plot), rects(old)
    assert [len(r) for r in new_r] == [adata.n_obs] * len(GENES)
    assert len(new_r) == len(old_r)
    for a, b in zip(new_r, old_r):
        assert np.allclose(a, b)
    if not use_raw:
        assert (data["value"] < 0).any()
        # a negative value hangs below the zero baseline
        assert any(r[:, 1].min() < 0 for r in new_r)


def test_P1_tracksplot_cell_order_follows_group(adata):
    d = ag.plot_tracksplot(adata, GENES[:1], GROUP, use_raw=True).data
    d = d.sort_values("cell_rank")
    codes = d[GROUP].cat.codes.to_numpy()
    assert (np.diff(codes) >= 0).all()


# --------------------------------------------------------------------------- #
# P4: ridge cardinality guard
# --------------------------------------------------------------------------- #
def test_P4_ridge_guard_names_field_and_level_count(adata):
    a = adata.copy()
    a.obs["continuous"] = pd.Categorical(np.arange(a.n_obs).astype(str))
    with pytest.raises(ValueError, match=r"'continuous' has 700 levels.*at most 100"):
        ag.plot_ridge(a, GENES[:1], "continuous", use_raw=True)


# --------------------------------------------------------------------------- #
# F2a/F2b: repeatable jitter
# --------------------------------------------------------------------------- #
SEEDED = {
    "box": lambda a, **k: ag.plot_box(a, GENES[:1], GROUP, use_raw=True, **k),
    "violin": lambda a, **k: ag.plot_violin(
        a, GENES[:1], GROUP, use_raw=True, add_points=True, **k
    ),
    "sina": lambda a, **k: ag.plot_sina(a, GENES[:1], GROUP, use_raw=True, **k),
}


@pytest.mark.parametrize("name", SEEDED)
@pytest.mark.parametrize("seed", [0, 3])
def test_F2_seeded_jitter_repeats_on_redraw_and_fresh_plots(adata, name, seed):
    plot = SEEDED[name](adata, random_state=seed)
    first, second = _png(plot), _png(plot)  # repeated exports of one plot
    assert first == second
    assert first == _png(SEEDED[name](adata, random_state=seed))  # fresh plot
    assert first != _png(SEEDED[name](adata, random_state=seed + 1))


POINT_LAYER = {"box": "geom_jitter", "violin": "geom_jitter", "sina": "geom_sina"}


def _point_x(plot, geom: str) -> np.ndarray:
    """Final x of the jittered layer after a draw (a draw rebuilds layer data)."""
    import matplotlib.pyplot as plt

    plt.close(plot.draw())
    layer = [lay for lay in plot.layers if type(lay.geom).__name__ == geom][0]
    return layer.data["x"].to_numpy().copy()


@pytest.mark.parametrize("name", SEEDED)
@pytest.mark.parametrize("seed", [0, 3])
def test_F2_point_positions_identical_across_draws(adata, name, seed):
    plot = SEEDED[name](adata, random_state=seed)
    geom = POINT_LAYER[name]
    first, second = _point_x(plot, geom), _point_x(plot, geom)  # state not consumed
    assert len(first) == adata.n_obs
    assert np.array_equal(first, second)
    assert np.array_equal(first, _point_x(SEEDED[name](adata, random_state=seed), geom))
    assert not np.array_equal(first, _point_x(SEEDED[name](adata, random_state=seed + 1), geom))


@pytest.mark.parametrize("name", SEEDED)
def test_F2_none_is_explicitly_nondeterministic(adata, name):
    a = _png(SEEDED[name](adata, random_state=None))
    b = _png(SEEDED[name](adata, random_state=None))
    assert a != b
    plot = SEEDED[name](adata, random_state=None)
    _png(plot)
    _png(plot)  # a second export of a None-seeded plot must not crash


# --------------------------------------------------------------------------- #
# F3a/F3b: points are visible and counted
# --------------------------------------------------------------------------- #
def test_F3a_violin_add_points_adds_one_mark_per_cell(adata):
    without = ag.plot_violin(adata, GENES[:1], GROUP, use_raw=True)
    with_pts = ag.plot_violin(adata, GENES[:1], GROUP, use_raw=True, add_points=True)
    assert sum(_scatter_counts(with_pts)) - sum(_scatter_counts(without)) == adata.n_obs
    jitter = [layer for layer in with_pts.layers if type(layer.geom).__name__ == "geom_jitter"]
    assert jitter[0].geom.aes_params["size"] >= 0.6 and jitter[0].geom.aes_params["alpha"] >= 0.4


def test_F3b_box_jitter_and_sina_points_visible(adata):
    box = ag.plot_box(adata, GENES[:1], GROUP, use_raw=True)
    jit = [layer for layer in box.layers if type(layer.geom).__name__ == "geom_jitter"][0]
    assert jit.geom.aes_params["size"] >= 0.6 and jit.geom.aes_params["alpha"] >= 0.4
    sina = ag.plot_sina(adata, GENES[:1], GROUP, use_raw=True)
    assert sum(_scatter_counts(sina)) == adata.n_obs
    sp = [layer for layer in sina.layers if type(layer.geom).__name__ == "geom_sina"][0]
    assert sp.geom.aes_params["size"] >= 0.8
    # legend keys are drawn larger and opaque than the data points
    assert sina.guides.color.override_aes == {"size": 3, "alpha": 1}


# --------------------------------------------------------------------------- #
# F5: statistical labels clear the violin tips / whiskers
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("name", ["violin", "box"])
def test_F5_stat_label_sits_above_the_data(adata, name):
    import matplotlib.pyplot as plt

    fn = {"violin": ag.plot_violin, "box": ag.plot_box}[name]
    plot = fn(adata, GENES, GROUP, use_raw=True, stats=True)
    fig = plot.draw()
    values = plot.data.groupby("feature", observed=True)["value"].max()
    for ax, top in zip([a for a in fig.axes if a.texts], values):
        label_y = max(t.get_position()[1] for t in ax.texts)
        assert label_y > top, "label overlaps the tallest violin tip / whisker"
        assert ax.get_ylim()[1] > label_y
    plt.close(fig)


# --------------------------------------------------------------------------- #
# F1 + U5: integer-data sina y-jitter is documented, never hidden by editing data
# --------------------------------------------------------------------------- #
def test_F1_U5_sina_does_not_alter_integer_values_and_documents_limit():
    import anndata as ad

    rng = np.random.RandomState(0)
    counts = rng.poisson(2, size=(120, 1)).astype(float)
    a = ad.AnnData(counts)
    a.var_names = ["g"]
    a.obs["grp"] = pd.Categorical(np.repeat(["a", "b", "c"], 40))
    plot = ag.plot_sina(a, ["g"], "grp")
    assert np.array_equal(plot.data["value"].to_numpy(), counts[:, 0])
    assert "integer" in ag.plot_sina.__doc__ and "U5" in ag.plot_sina.__doc__


# --------------------------------------------------------------------------- #
# Q1: no pandas categories-order deprecation
# --------------------------------------------------------------------------- #
def test_Q1_categories_order_emits_no_deprecation(adata):
    cats = list(adata.obs[GROUP].cat.categories)[::-1]
    calls = [
        lambda: ag.plot_violin(adata, GENES, GROUP, categories_order=cats),
        lambda: ag.plot_box(adata, GENES, GROUP, categories_order=cats),
        lambda: ag.plot_sina(adata, GENES, GROUP, categories_order=cats),
        lambda: ag.plot_ridge(adata, GENES, GROUP, categories_order=cats, n_grid=16),
        lambda: ag.plot_expression_bar(adata, GENES, GROUP, categories_order=cats),
        lambda: ag.plot_stacked_violin(adata, GENES, GROUP, categories_order=cats),
        lambda: ag.plot_proportions(adata, GROUP, "depth", categories_order=cats),
    ]
    with warnings.catch_warnings():
        warnings.filterwarnings("error", message=".*[Cc]ategor.*")
        for call in calls:
            call().draw()


# --------------------------------------------------------------------------- #
# Q3, Q13: density
# --------------------------------------------------------------------------- #
def test_Q3_one_feature_joint_warns(adata):
    pytest.importorskip("pynebulosa")
    with pytest.warns(UserWarning, match="joint=True needs at least two features"):
        plot = ag.plot_density(adata, "CD3D", joint=True, use_raw=True)
    assert list(plot.data["feature"].cat.categories) == ["CD3D"]
    assert "absolute density" in ag.plot_density.__doc__


def test_Q13_density_downsample_and_seed(adata):
    pytest.importorskip("pynebulosa")
    a = ag.plot_density(adata, "CD3D", use_raw=True, downsample=100, random_state=1)
    b = ag.plot_density(adata, "CD3D", use_raw=True, downsample=100, random_state=1)
    c = ag.plot_density(adata, "CD3D", use_raw=True, downsample=100, random_state=2)
    assert len(a.data) == 100
    assert a.data.equals(b.data)
    assert not a.data.equals(c.data)
    with pytest.raises(ValueError, match="downsample"):
        ag.plot_density(adata, "CD3D", use_raw=True, downsample=0)


# --------------------------------------------------------------------------- #
# Q4: sina bins and fill-legend warning
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("bad", [0, -1, 2.5, True])
def test_Q4_sina_bins_validated(adata, bad):
    with pytest.raises(ValueError, match="bins"):
        ag.plot_sina(adata, GENES[:1], GROUP, bins=bad)


def test_Q4_sina_without_violin_has_no_fill_warning(adata):
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        ag.plot_sina(adata, GENES[:1], GROUP, use_raw=True, violin=False).draw()
    assert not [
        w
        for w in caught
        if issubclass(w.category, PlotnineWarning) and "fill" in str(w.message).lower()
    ]


# --------------------------------------------------------------------------- #
# Q12: ridge blank panels and bandwidth
# --------------------------------------------------------------------------- #
def test_Q12_ridge_warns_for_constant_gene_and_bandwidth_changes_curve(adata):
    a = adata.copy()
    a.obs["flat"] = 0.0
    with pytest.warns(UserWarning, match=r"constant values for \['flat'\]"):
        ag.plot_ridge(a, ["flat"], GROUP).draw()
    narrow = ag.plot_ridge(adata, GENES[:1], GROUP, use_raw=True, bw_adjust=0.3).data
    wide = ag.plot_ridge(adata, GENES[:1], GROUP, use_raw=True, bw_adjust=3.0).data
    assert not np.allclose(narrow["ymax"], wide["ymax"])
    with pytest.raises(ValueError, match="bw_adjust"):
        ag.plot_ridge(adata, GENES[:1], GROUP, bw_adjust=0)


def test_Q12_stacked_violin_hides_crowded_y_labels_by_default(adata):
    genes = list(adata.raw.var_names[:10])
    for swap in (False, True):
        plot = ag.plot_stacked_violin(adata, genes, GROUP, use_raw=True, swap_axes=swap)
        fig = plot.draw()
        assert all(not [t for t in ax.get_yticklabels() if t.get_text()] for ax in fig.axes)
        kw = dict(use_raw=True, yticklabels=True, swap_axes=swap)
        shown = ag.plot_stacked_violin(adata, genes[:2], GROUP, **kw).draw()
        assert any(t.get_text() for ax in shown.axes for t in ax.get_yticklabels())


# --------------------------------------------------------------------------- #
# Q16: expression-line bins, split-violin ncol, identical group/split
# --------------------------------------------------------------------------- #
def test_Q16_expression_line_bins_numeric_x(adata):
    plain = ag.plot_expression_line(adata, GENES[:1], "n_genes", use_raw=True)
    binned = ag.plot_expression_line(adata, GENES[:1], "n_genes", use_raw=True, bins=5)
    assert binned.data["n_genes"].nunique() == 5
    assert plain.data["n_genes"].nunique() > 5
    lo, hi = adata.obs["n_genes"].min(), adata.obs["n_genes"].max()
    assert lo <= binned.data["n_genes"].min() and binned.data["n_genes"].max() <= hi
    binned.draw()
    with pytest.raises(ValueError, match="bins"):
        ag.plot_expression_line(adata, GENES[:1], "n_genes", bins=0)
    with pytest.raises(ValueError, match="numeric"):
        ag.plot_expression_line(adata, GENES[:1], "depth", bins=3)


def test_Q16_split_violin_ncol_is_ignored_with_split_by(adata):
    plot = ag.plot_violin(adata, GENES, GROUP, split_by="depth", ncol=2, use_raw=True)
    assert type(plot.facet).__name__ == "facet_grid"
    assert "ncol" in ag.plot_violin.__doc__ and "ignored" in ag.plot_violin.__doc__
    plot.draw()


def test_Q16_proportions_identical_group_and_split_rejected(adata):
    with pytest.raises(ValueError, match="both 'bulk_labels'"):
        ag.plot_proportions(adata, GROUP, GROUP)


# --------------------------------------------------------------------------- #
# Q9: tracksplot group labels on x
# --------------------------------------------------------------------------- #
def test_Q9_tracksplot_labels_each_group_at_block_centre(adata):
    plot = ag.plot_tracksplot(adata, GENES, GROUP, use_raw=True)
    fig = plot.draw()
    labels = [t.get_text() for t in fig.axes[-1].get_xticklabels()]
    assert labels == list(plot.data[GROUP].cat.categories)


# --------------------------------------------------------------------------- #
# VERIFY-4: non-mutation, later `+`, tracksplot visible population
# --------------------------------------------------------------------------- #
HELPERS = {
    "violin": lambda a: ag.plot_violin(a, GENES, GROUP, use_raw=True, add_points=True, stats=True),
    "box": lambda a: ag.plot_box(a, GENES, GROUP, use_raw=True, stats=True),
    "sina": lambda a: ag.plot_sina(a, GENES, GROUP, use_raw=True, downsample=50),
    "stacked": lambda a: ag.plot_stacked_violin(a, GENES, GROUP, use_raw=True, stripplot=True),
    "ridge": lambda a: ag.plot_ridge(a, GENES, GROUP, use_raw=True, n_grid=16),
    "bar": lambda a: ag.plot_expression_bar(a, GENES, GROUP, use_raw=True),
    "line": lambda a: ag.plot_expression_line(
        a, GENES, "n_genes", group_by=GROUP, use_raw=True, bins=4
    ),
    "tracksplot": lambda a: ag.plot_tracksplot(a, GENES, GROUP, use_raw=True),
    "obs_violin": lambda a: ag.plot_violin(a, ["n_genes", ag.obs("n_counts")], GROUP, log=True),
}


@pytest.mark.parametrize("name", HELPERS)
def test_VERIFY4_distribution_helpers_do_not_mutate_adata(adata, name):
    a = adata.copy()
    before = _key(a)
    HELPERS[name](a).draw()
    assert _key(a) == before


@pytest.mark.parametrize("name", HELPERS)
def test_VERIFY4_later_plus_modifications_still_build(adata, name):
    plot = HELPERS[name](adata) + labs(title="later") + theme(legend_position="bottom")
    if name not in ("tracksplot", "stacked"):
        plot = plot + scale_y_continuous(limits=(None, None))
    plot.draw()
    assert plot.labels.title == "later"


def test_VERIFY4_tracksplot_missing_group_not_visible(adata):
    import matplotlib.pyplot as plt

    a = adata.copy()
    a.obs[GROUP] = a.obs[GROUP].astype(object)
    a.obs.loc[a.obs_names[:25], GROUP] = np.nan
    a.obs[GROUP] = pd.Categorical(a.obs[GROUP])
    with pytest.warns(UserWarning, match="dropped 25 of 700"):
        plot = ag.plot_tracksplot(a, GENES, GROUP, use_raw=True)
    fig = plot.draw()
    n_rects = [len(c.get_paths()) for ax in fig.axes for c in ax.collections]
    assert n_rects == [a.n_obs - 25] * len(GENES)
    legend_text = " ".join(t.get_text() for t in fig.findobj(lambda o: hasattr(o, "get_text")))
    assert "nan" not in legend_text.lower().split()
    plt.close(fig)


@pytest.mark.parametrize("swap", [False, True])
def test_stacked_violin_first_x_label_inside_canvas(adata, swap):
    # hidden y labels must not take away the space the first rotated x label needs
    from matplotlib.backends.backend_agg import FigureCanvasAgg

    genes = ["CD3D", "NKG7", "CST3"]
    plot = ag.plot_stacked_violin(adata, genes, "bulk_labels", use_raw=True, swap_axes=swap)
    fig = plot.draw()  # default rotation: 45 for genes, 90 for groups
    canvas = FigureCanvasAgg(fig)
    canvas.draw()
    renderer = canvas.get_renderer()
    left = min(
        t.get_window_extent(renderer).x0
        for ax in fig.axes
        for t in ax.get_xticklabels()
        if t.get_visible() and t.get_text()
    )
    assert left >= 0
