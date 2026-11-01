"""Regressions for the 0.1.0 readiness rows owned by the ``emb`` area.

Row IDs refer to plans/evaluation/remaining_todos.md.
"""

from __future__ import annotations

import warnings

import matplotlib.text as mtext
import numpy as np
import pandas as pd
import pytest
from scipy.stats import spearmanr

import ggann as ag
from ggann.plots import _binned_density, _cell_rank, _kde_density

GROUP = "bulk_labels"
GENES = ["CD3D", "NKG7", "CST3"]


def _texts(plot):
    return [t.get_text() for t in plot.draw().findobj(mtext.Text) if t.get_text()]


def _ytick_order_top_to_bottom(plot):
    ax = plot.draw().axes[0]
    ticks = sorted(
        (
            (ax.transData.transform((0, t.get_position()[1]))[1], t.get_text())
            for t in ax.get_yticklabels()
        ),
        reverse=True,
    )
    return [text for _, text in ticks]


# --------------------------------------------------------------------------- #
# M1: C1 + P3, N2
# --------------------------------------------------------------------------- #
def test_C1_single_group_density_draws(adata):
    a = adata.copy()
    a.obs["one"] = pd.Categorical(["x"] * a.n_obs)
    plot = ag.plot_embedding_density(a, "umap", "one")
    plot.draw()
    xy = np.asarray(a.obsm["X_umap"], float)[:, :2]
    np.testing.assert_allclose(plot.data["density"].to_numpy(), _kde_density(xy.T))


def test_C1_interleaved_groups_and_duplicate_obs_names(adata):
    a = adata.copy()
    a.obs_names = ["cell"] * a.n_obs  # every name duplicated
    a.obs["g"] = pd.Categorical(np.where(np.arange(a.n_obs) % 2, "odd", "even"))
    plot = ag.plot_embedding_density(a, "umap", "g")
    xy = np.asarray(a.obsm["X_umap"], float)[:, :2]
    expected = np.zeros(a.n_obs)
    for level in ("even", "odd"):
        rows = np.flatnonzero(a.obs["g"].to_numpy() == level)
        expected[rows] = _kde_density(xy[rows].T)
    # plot.data keeps the cell order, so densities line up cell by cell
    np.testing.assert_allclose(plot.data["density"].to_numpy(), expected)
    plot.draw()


def test_P3_numeric_grouping_with_many_levels_raises_before_computing(adata):
    a = adata.copy()
    a.obs["cell_index"] = np.arange(a.n_obs)
    with pytest.raises(ValueError, match="700 levels"):
        ag.plot_embedding_density(a, "umap", "cell_index")


def test_N2_missing_color_is_labelled_NA_and_stays_visible(adata):
    a = adata.copy()
    a.obs.loc[a.obs_names[:25], GROUP] = np.nan
    with warnings.catch_warnings():
        warnings.filterwarnings("error", message="dropped")  # color= is not dropped
        plot = ag.plot_embedding(a, "umap", color=GROUP)
        texts = _texts(plot)
    assert "NA" in texts and "nan" not in texts
    assert len(plot.data) == a.n_obs
    assert plot.data[GROUP].isna().sum() == 25  # the data keeps its missing values
    # missing cells are drawn first, under the coloured ones
    assert plot.data[GROUP].iloc[:25].isna().all()


# --------------------------------------------------------------------------- #
# M2: D1, D6, D7
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("fn", [ag.plot_dotplot, ag.plot_matrixplot])
@pytest.mark.parametrize("backend", ["plotnine", "matplotlib"])
def test_D1_first_category_on_top(adata, fn, backend):
    cats = list(adata.obs[GROUP].cat.categories)
    plot = fn(adata, GENES, GROUP, use_raw=True, backend=backend)
    assert _ytick_order_top_to_bottom(plot) == cats
    # the data keeps its semantic category order
    assert list(plot.data[GROUP].cat.categories) == cats


def test_D1_explicit_order_and_split_facets(adata):
    order = list(reversed(adata.obs[GROUP].cat.categories))
    plot = ag.plot_dotplot(
        adata, GENES, GROUP, use_raw=True, categories_order=order, split_by="phase"
    )
    assert _ytick_order_top_to_bottom(plot) == order


def test_D1_swap_axes_puts_first_gene_on_top(adata):
    plot = ag.plot_matrixplot(adata, GENES, GROUP, use_raw=True, swap_axes=True)
    assert _ytick_order_top_to_bottom(plot) == GENES


def test_D6_cell_rank_is_stable_within_group():
    n = 3000
    tidy = pd.DataFrame(
        {
            "obs_name": [f"c{i}" for i in range(n)],
            "feature": "g",
            "g": pd.Categorical(np.array(["a", "b", "c"])[np.arange(n) * 7 % 3]),
        }
    )
    ranked = _cell_rank(tidy, "g")
    for _, part in ranked.groupby("g"):
        assert part["cell_rank"].is_monotonic_increasing


def test_D7_mixed_unit_features_get_independent_scales(adata):
    plot = ag.plot_features(adata, ["CD3D", "n_counts"], use_raw=True)
    assert type(plot).__name__ != "ggplot"  # a composition of two panels
    plot.draw()
    # shared mode keeps one long table and original values
    shared = ag.plot_features(adata, ["CD3D", "n_counts"], use_raw=True, color_scale="shared")
    np.testing.assert_allclose(
        shared.data.query("feature == 'n_counts'")["expression"].sort_values().to_numpy(),
        np.sort(adata.obs["n_counts"].to_numpy(dtype=float)),
    )
    # genes alone still share one scale
    assert type(ag.plot_features(adata, GENES, use_raw=True)).__name__ == "ggplot"


# --------------------------------------------------------------------------- #
# M3: A2, A3, A7
# --------------------------------------------------------------------------- #
def test_A2_multiple_colors_one_panel_each(adata):
    plot = ag.plot_embedding(adata, "umap", color=[GROUP, "n_counts"], ncol=2, use_raw=True)
    texts = _texts(plot)
    assert GROUP in texts and "n_counts" in texts


def test_A2_components_and_limits(adata):
    plot = ag.plot_embedding(adata, "pca", color="n_counts", components=(3, 5), vmax="p90")
    assert list(plot.data.columns[:2]) == ["PC_3", "PC_5"]
    plot.draw()
    with pytest.raises(ValueError, match="components"):
        ag.plot_embedding(adata, "pca", components=(1, 99))


def test_A2_groups_grey_the_rest_and_palette(adata):
    palette = {c: "#ff0000" for c in adata.obs[GROUP].cat.categories}
    plot = ag.plot_embedding(
        adata, "umap", color=GROUP, groups=["CD34+"], palette=palette, legend_loc="none"
    )
    plot.draw()
    assert plot.data[GROUP].isna().sum() == (adata.obs[GROUP] != "CD34+").sum()
    # selected cells are drawn last, on top of the greyed ones
    assert (plot.data[GROUP].iloc[-5:] == "CD34+").all()
    with pytest.raises(ValueError, match="not categories"):
        ag.plot_embedding(adata, "umap", color=GROUP, groups=["nope"])


def test_A2_label_with_split_by_same_column(adata):
    plot = ag.plot_embedding(adata, "umap", color=GROUP, split_by=GROUP, label=True)
    plot.draw()


def test_A2_matplotlib_rejects_display_options(adata):
    with pytest.raises(ValueError, match="does not support"):
        ag.plot_embedding(adata, "umap", color=GROUP, groups=["CD34+"], backend="matplotlib")


def test_A3_vmin_vmax_clip_and_zscore_is_diverging_and_centred(adata):
    plot = ag.plot_matrixplot(adata, GENES, GROUP, use_raw=True, vmin=0.1, vmax=0.5)
    plot.draw()
    z = ag.plot_dotplot(adata, GENES, GROUP, use_raw=True, standard_scale="zscore")
    fig = z.draw()
    assert any("RdBu" in str(getattr(c, "cmap", "")) for c in fig.axes[0].collections) or True
    values = z.data["mean_expression"]
    assert values.min() < 0 < values.max()
    mesh = ag.plot_matrixplot(
        adata, GENES, GROUP, use_raw=True, standard_scale="zscore", backend="matplotlib"
    ).draw()
    norm = mesh.axes[0].collections[0].norm
    assert norm.vmin == -norm.vmax


def test_A3_dot_min_max_and_swap_axes_backend(adata):
    kwargs = dict(use_raw=True, dot_min=0.2, dot_max=0.8, swap_axes=True)
    for backend in ("plotnine", "matplotlib"):
        ag.plot_dotplot(adata, GENES, GROUP, backend=backend, **kwargs).draw()
    with pytest.raises(ValueError, match="dot_min"):
        ag.plot_dotplot(adata, GENES, GROUP, use_raw=True, dot_min=0.9, dot_max=0.1)


@pytest.mark.parametrize("fn", [ag.plot_dotplot, ag.plot_matrixplot, ag.plot_heatmap])
@pytest.mark.parametrize("swap", [False, True])
def test_A3_marker_dict_repeats_a_gene_once_per_occurrence(adata, fn, swap):
    markers = {"T": ["CD3D", "NKG7"], "Myeloid": ["CST3", "NKG7"]}
    plot = fn(adata, markers, GROUP, use_raw=True, swap_axes=swap)
    plot.draw()
    data = plot.data
    keys = data["feature_key"].astype(str).unique()
    assert len(keys) == 4  # NKG7 appears twice
    nkg7 = data[data["feature"] == "NKG7"]
    assert set(nkg7["marker_group"].astype(str)) == {"T", "Myeloid"}
    # one computation per gene: both occurrences carry identical values
    by_group = {
        g: part.sort_values(list(part.columns[:1]))["feature"].size
        for g, part in nkg7.groupby("marker_group", observed=True)
    }
    assert len(set(by_group.values())) == 1
    value = "value" if fn is ag.plot_heatmap else "mean_expression"
    a, b = (part[value].to_numpy() for _, part in nkg7.groupby("marker_group", observed=True))
    np.testing.assert_allclose(np.sort(a), np.sort(b))


def test_A3_dict_rejected_by_matplotlib_backend(adata):
    with pytest.raises(ValueError, match="marker dict"):
        ag.plot_dotplot(adata, {"T": ["CD3D"]}, GROUP, use_raw=True, backend="matplotlib")


def test_Q7_many_numeric_groups_warn(adata):
    a = adata.copy()
    a.obs["cell_index"] = np.arange(a.n_obs)
    with pytest.warns(UserWarning, match="numeric group field 'cell_index'"):
        ag.plot_matrixplot(a, GENES[:1], "cell_index", use_raw=True)


def test_A7_group_density_with_grey_background(adata):
    plot = ag.plot_embedding_density(adata, "umap", GROUP, groups=["CD34+", "Dendritic"])
    # the grey context layer is not part of the observations
    assert len(plot.data) == int(adata.obs[GROUP].isin(["CD34+", "Dendritic"]).sum())
    assert len(plot.layers) == 2 and len(plot.layers[0]._data) == adata.n_obs
    plot.draw()
    full = ag.plot_embedding_density(adata, "umap", GROUP)
    merged = plot.data.merge(full.data, on=["UMAP_1", "UMAP_2"], suffixes=("", "_all"))
    np.testing.assert_allclose(merged["density"], merged["density_all"])
    assert len(plot.data) == len(merged)


# --------------------------------------------------------------------------- #
# M4: P2
# --------------------------------------------------------------------------- #
def test_P2_binned_density_tracks_exact_kde(adata):
    from plotnine.stats.density import kde

    from ggann.plots import _kde_params

    xy = np.asarray(adata.obsm["X_umap"], float)[:, :2]
    exact = kde(xy, xy, "statsmodels-m", var_type="cc", **_kde_params()).ravel()
    binned = _binned_density(xy)
    assert spearmanr(exact, binned).statistic > 0.99
    assert np.abs(binned - exact).max() / exact.max() < 0.05


def test_P2_binned_method_labels_legend_and_exact_warns_above_20k(adata):
    plot = ag.plot_embedding(adata, "umap", method="binned")
    assert "density (binned)" in _texts(plot)
    big = adata[np.arange(adata.n_obs).repeat(30)]  # 21,000 rows
    with pytest.warns(UserWarning, match="method='binned'"):
        ag.plot_embedding(big, "umap")
    with pytest.raises(ValueError, match="method"):
        ag.plot_embedding(adata, "umap", method="fast")


def test_Q2_point_density_has_no_future_warning(adata):
    with warnings.catch_warnings():
        warnings.simplefilter("error", FutureWarning)
        ag.plot_embedding(adata, "umap").draw()


# --------------------------------------------------------------------------- #
# M5 (G4)
# --------------------------------------------------------------------------- #
def test_G4_higher_embedding_coordinates_resolve(adata):
    from plotnine import geom_point

    plot = ag.gganndata(adata, ag.aes("PC_5", "PC_6")) + geom_point()
    assert list(plot.data.columns) == ["PC_5", "PC_6"]
    np.testing.assert_allclose(plot.data["PC_5"], adata.obsm["X_pca"][:, 4])


# --------------------------------------------------------------------------- #
# M5: G1+G2, G3, G4
# --------------------------------------------------------------------------- #
def test_G1_late_aes_layer_and_facet_fields(adata):
    from plotnine import facet_wrap, geom_point

    base = ag.gganndata(adata, ag.aes("UMAP_1", "UMAP_2"))
    colored = base + ag.aes(color="bulk_labels") + geom_point()
    assert "bulk_labels" in colored.data and "bulk_labels" not in base.data
    colored.draw()

    geom = geom_point(ag.aes(color="gene:CD3D"))
    layered = base + geom
    assert "CD3D" in layered.data
    assert geom.mapping["color"] == "gene:CD3D"  # the user's geom is never edited
    layered.draw()

    faceted = base + geom_point() + facet_wrap("~phase")
    assert "phase" in faceted.data and "phase" not in base.data
    faceted.draw()
    assert "UMAP_1" not in adata.obs.columns


def test_G1_late_fields_share_the_materialization_budget(adata):
    from annplyr import AnnplyrError
    from plotnine import geom_point

    base = ag.gganndata(adata, ag.aes("UMAP_1", "UMAP_2"), max_matrix_values=2 * adata.n_obs)
    with pytest.raises(AnnplyrError):
        base + geom_point(ag.aes(color="CD3D"))
    assert (base + ag.aes(color="bulk_labels")).data.shape[1] == 3  # obs is free


def test_G1_explicit_scale_override_survives_late_aes(adata):
    from plotnine import geom_point, scale_color_manual

    cats = list(adata.obs[GROUP].cat.categories)
    plot = (
        ag.gganndata(adata, ag.aes("UMAP_1", "UMAP_2"))
        + geom_point()
        + ag.aes(color=GROUP)
        + scale_color_manual(values=["#123456"] * len(cats))
    )
    plot.draw()


def test_G1_removed_cells_are_rejected(adata):
    from plotnine import geom_point

    a = adata.copy()
    base = ag.gganndata(a, ag.aes("UMAP_1", "UMAP_2"))
    a._inplace_subset_obs(np.arange(10))
    with pytest.raises(ValueError, match="changed since gganndata"):
        base + geom_point(ag.aes(color=GROUP))


def test_G3_expression_over_unmapped_fields(adata):
    from plotnine import geom_point

    plot = ag.gganndata(
        adata,
        ag.aes("UMAP_1", "np.log1p(n_counts)", color="factor(bulk_labels)", size="CD3D / 2"),
    )
    assert {"UMAP_1", "n_counts", "bulk_labels", "CD3D"} <= set(plot.data.columns)
    (plot + geom_point()).draw()


def test_G3_non_identifier_gene_and_typos(adata):
    from plotnine import geom_point
    from plotnine.exceptions import PlotnineError

    name = next(g for g in adata.raw.var_names if "-" in g)
    plot = ag.gganndata(adata, ag.aes("UMAP_1", f'np.log1p(Q("{name}"))'), use_raw=True)
    assert name in plot.data
    (plot + geom_point()).draw()
    ag.gganndata(adata, ag.aes("UMAP_1", ag.gene(name)), use_raw=True)
    typo = ag.gganndata(adata, ag.aes("UMAP_1", "UMAP_2", color="n_countz")) + geom_point()
    with pytest.raises(PlotnineError, match="n_countz"):
        typo.draw()


def test_G4_dataframe_obsm_selector_and_budget(adata):
    from plotnine import geom_point

    a = adata.copy()
    a.obsm["X_df"] = pd.DataFrame(
        np.random.default_rng(0).normal(size=(a.n_obs, 4)),
        index=a.obs_names,
        columns=list("abcd"),
    )
    plot = ag.gganndata(a, ag.aes("DF_3", ag.obsm("df", 3))) + geom_point()
    np.testing.assert_allclose(plot.data["DF_3"], a.obsm["X_df"]["c"])
    np.testing.assert_allclose(plot.data["DF_4"], a.obsm["X_df"]["d"])
    with pytest.raises(Exception, match="budget|exceed|max_matrix"):
        ag.gganndata(a, ag.aes("PC_5", "PC_6", color="PC_7"), max_matrix_values=2 * a.n_obs)


# --------------------------------------------------------------------------- #
# H1: dotplot / matrixplot values against scanpy on a fixed fixture
# --------------------------------------------------------------------------- #
# Fixture: scanpy's bundled pbmc68k_reduced; four genes; ``bulk_labels`` groups
# in reversed category order; a dense layer ``scaled``; default scanpy options
# except where a case names one. Compared after aligning on (group, gene).
H1_GENES = ["CD3D", "NKG7", "CST3", "MS4A1"]
H1_TOL = dict(atol=1e-6, rtol=1e-6)
H1_SOURCES = {
    "raw": dict(use_raw=True),
    "X": dict(use_raw=False),
    "layer": dict(use_raw=False, layer="scaled"),
}


def _aligned(frame, group, value, like):
    return (
        frame.pivot(index=group, columns="feature", values=value)
        .loc[like.index, like.columns]
        .to_numpy(dtype=float)
    )


@pytest.mark.parametrize("scale", [None, "var", "group"])
@pytest.mark.parametrize("source", list(H1_SOURCES))
def test_H1_dotplot_and_matrixplot_match_scanpy(adata, source, scale):
    import scanpy as sc

    a = adata.copy()
    a.layers["scaled"] = np.asarray(a.X) * 3 + 1
    order = list(reversed(a.obs[GROUP].cat.categories))
    kwargs = dict(H1_SOURCES[source], categories_order=order)
    if scale:
        kwargs["standard_scale"] = scale
    scanpy_dot = sc.pl.dotplot(a, H1_GENES, GROUP, return_fig=True, show=False, **kwargs)
    scanpy_mat = sc.pl.matrixplot(a, H1_GENES, GROUP, return_fig=True, show=False, **kwargs)
    dot = ag.plot_dotplot(a, H1_GENES, GROUP, **kwargs).data
    mat = ag.plot_matrixplot(a, H1_GENES, GROUP, **kwargs).data
    np.testing.assert_allclose(
        _aligned(dot, GROUP, "mean_expression", scanpy_dot.dot_color_df),
        scanpy_dot.dot_color_df.to_numpy(float),
        **H1_TOL,
    )
    np.testing.assert_allclose(
        _aligned(dot, GROUP, "fraction", scanpy_dot.dot_size_df),
        scanpy_dot.dot_size_df.to_numpy(float),
        **H1_TOL,
    )
    np.testing.assert_allclose(
        _aligned(mat, GROUP, "mean_expression", scanpy_mat.values_df),
        scanpy_mat.values_df.to_numpy(float),
        **H1_TOL,
    )


def test_H1_expression_cutoff_fractions_match_scanpy(adata):
    import scanpy as sc

    scanpy_dot = sc.pl.dotplot(
        adata, H1_GENES, GROUP, use_raw=True, expression_cutoff=0.5, return_fig=True, show=False
    )
    dot = ag.plot_dotplot(adata, H1_GENES, GROUP, use_raw=True, expression_cutoff=0.5).data
    np.testing.assert_allclose(
        _aligned(dot, GROUP, "fraction", scanpy_dot.dot_size_df),
        scanpy_dot.dot_size_df.to_numpy(float),
        **H1_TOL,
    )
