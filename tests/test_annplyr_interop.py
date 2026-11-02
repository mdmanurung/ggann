"""annplyr pipelines feed ggann without adapters.

Each test builds an annplyr result (a ``TidyAnnData`` or a table) and hands it
straight to ggann or plotnine, then checks the plotted data against the same
numbers read from the original AnnData. Uses the bundled ``pbmc68k_reduced``
through the shared ``adata`` fixture.
"""

from __future__ import annotations

import copy

import annplyr as ap
import numpy as np
import pandas as pd
import pytest
from annplyr.data import X, X_umap, obs, raw
from plotnine import aes, geom_boxplot, geom_col, geom_point, ggplot

import ggann as ag


def _layer_data(plot, i=0):
    # plotnine<0.16 has no ggplot.layer_data; build a copy and read the layer
    built = copy.deepcopy(plot)
    built._build()
    return built.layers[i].data


@pytest.fixture(scope="module")
def tidy(adata):
    return ap.tidy(adata)


def raw_gene(adata, gene):
    column = adata.raw[:, gene].X
    return np.asarray(column.toarray() if hasattr(column, "toarray") else column).ravel()


# ---- TidyAnnData results into gganndata and plot_* helpers ------------------


def test_filter_result_feeds_gganndata(adata, tidy):
    kept = tidy.filter(obs.n_counts > 3333)
    plot = ag.gganndata(kept, aes("UMAP_1", "UMAP_2", color="louvain")) + geom_point()
    assert 0 < kept.n_obs < adata.n_obs
    assert list(plot.data.index) == list(kept.obs_names)
    np.testing.assert_allclose(plot.data["UMAP_1"], kept.obsm["X_umap"][:, 0])
    assert list(plot.data["louvain"].astype(str)) == list(kept.obs["louvain"].astype(str))


def test_mutate_column_in_aes_and_color(adata, tidy):
    scored = tidy.mutate(b_score=(raw.MS4A1 + raw.CD79A) / 2)
    expected = (raw_gene(adata, "MS4A1") + raw_gene(adata, "CD79A")) / 2
    plot = ag.gganndata(scored, aes("UMAP_1", "UMAP_2", color="b_score")) + geom_point()
    np.testing.assert_allclose(plot.data["b_score"].astype(float), expected, rtol=1e-5)
    built = ag.plot_embedding(scored, color="b_score")
    # plot_embedding sorts cells by value so high scores draw on top: align by name
    aligned = built.data["b_score"].astype(float).loc[scored.obs_names]
    np.testing.assert_allclose(aligned, expected, rtol=1e-5)
    # a mutate column also works on the y axis of a grammar plot
    by_group = ag.gganndata(scored, aes("louvain", "b_score")) + geom_boxplot()
    by_group._build()
    assert _layer_data(by_group)["ymax"].max() <= expected.max() + 1e-6


def test_rename_result_is_a_plain_obs_column(tidy):
    renamed = tidy.rename(cluster=obs.louvain)
    plot = ag.plot_embedding(renamed, color="cluster")
    assert "cluster" in plot.data.columns
    assert list(plot.data["cluster"].astype(str)) == list(renamed.obs["cluster"].astype(str))


def test_slice_max_by_group_feeds_embedding(tidy):
    top = tidy.slice_max(obs.n_counts, n=3, by=obs.louvain)
    n_groups = tidy.obs["louvain"].nunique()
    assert top.n_obs == 3 * n_groups
    plot = ag.plot_embedding(top, color="louvain")
    assert len(plot.data) == top.n_obs
    assert plot.data["louvain"].value_counts().eq(3).all()


def test_aggregate_cells_pseudobulk_feeds_grammar_and_helpers(tidy):
    pseudobulk = tidy.aggregate_cells(obs.phase, obs.louvain, source=X)
    assert "n_cells" in pseudobulk.obs
    plot = ag.gganndata(pseudobulk, aes("louvain", "n_cells", fill="phase")) + geom_col()
    plot._build()
    assert plot.data["n_cells"].sum() == tidy.n_obs
    dot = ag.plot_dotplot(pseudobulk, ["CST3", "NKG7"], group_by="louvain")
    assert set(dot.data.columns) >= {"louvain"}
    assert len(dot.data) == 2 * pseudobulk.obs["louvain"].nunique()


def test_pipe_chain_returns_a_plot(tidy):
    plot = tidy.filter(obs.n_counts > 3333).pipe(ag.plot_embedding, color="louvain")
    assert isinstance(plot, ggplot)
    plot._build()
    assert _layer_data(plot)["x"].notna().all()


@pytest.mark.parametrize(
    "helper, args",
    [
        (ag.plot_violin, (["CST3"],)),
        (ag.plot_dotplot, (["CST3", "NKG7"],)),
        (ag.plot_matrixplot, (["CST3", "NKG7"],)),
    ],
)
def test_filtered_result_feeds_expression_helpers(tidy, helper, args):
    kept = tidy.filter(obs.n_counts > 3333)
    plot = helper(kept, *args, group_by="louvain")
    plot._build()
    assert set(plot.data["louvain"].astype(str)) <= set(kept.obs["louvain"].astype(str))


# ---- extraction tables into plain plotnine ----------------------------------


def test_to_wide_umap_names_draw_with_plain_plotnine(tidy):
    wide = tidy.to_wide(obs.louvain, X_umap)
    assert list(wide.columns) == ["louvain", "UMAP_1", "UMAP_2"]
    plot = ggplot(wide, aes("UMAP_1", "UMAP_2", color="louvain")) + geom_point() + ag.theme_ggann()
    plot._build()
    np.testing.assert_allclose(_layer_data(plot)["x"], tidy.obsm["X_umap"][:, 0])


def test_to_wide_names_match_the_embedding_ggann_resolves(tidy):
    wide = tidy.to_wide(X_umap)
    plot = ag.gganndata(tidy, aes("UMAP_1", "UMAP_2")) + geom_point()
    np.testing.assert_allclose(wide["UMAP_1"], plot.data["UMAP_1"])
    np.testing.assert_allclose(wide["UMAP_2"], plot.data["UMAP_2"])


def test_join_features_long_boxplot_needs_dense_abundance(adata, tidy):
    long = tidy.join_features(raw.CST3, raw.NKG7, shape="long")
    assert {".cell", ".feature", ".abundance_raw"} <= set(long.columns)
    assert isinstance(long[".abundance_raw"].dtype, pd.SparseDtype)
    sparse_plot = ggplot(long, aes("louvain", ".abundance_raw", fill=".feature")) + geom_boxplot()
    # plotnine cannot compute box statistics on a pandas sparse column
    with pytest.raises(ValueError, match="null value"):
        sparse_plot._build()
    dense = long.assign(**{".abundance_raw": long[".abundance_raw"].sparse.to_dense()})
    plot = ggplot(dense, aes("louvain", ".abundance_raw", fill=".feature")) + geom_boxplot()
    plot._build()
    cst3 = dense.loc[dense[".feature"] == "CST3", ".abundance_raw"].to_numpy()
    np.testing.assert_allclose(cst3, raw_gene(adata, "CST3"), rtol=1e-6)
    assert _layer_data(plot)["ymax"].max() <= dense[".abundance_raw"].max() + 1e-6


def test_wide_join_features_puts_genes_in_obs_for_grammar(adata, tidy):
    with_gene = tidy.join_features(raw.CST3, shape="wide")
    assert "CST3" in with_gene.obs.columns
    # the gene now lives in obs; ggann prefers obs and says so
    with pytest.warns(UserWarning, match="both an obs column and a gene"):
        plot = ag.gganndata(with_gene, aes("louvain", "CST3")) + geom_boxplot()
    plot._build()
    np.testing.assert_allclose(plot.data["CST3"].astype(float), raw_gene(adata, "CST3"), rtol=1e-6)


# ---- palette and category order ---------------------------------------------


def test_store_palette_colours_reach_the_drawn_points(tidy):
    levels = list(tidy.obs["louvain"].cat.categories)
    palette = {level: f"#{i * 20:02x}32c8" for i, level in enumerate(levels)}
    painted = ap.store_palette(tidy, "louvain", palette)
    assert ap.get_palette(painted, "louvain") == palette
    assert ag.obs_colors(painted, "louvain") == palette
    plot = ag.plot_embedding(painted, color="louvain")
    plot._build()
    drawn = _layer_data(plot)["color"].to_numpy()
    expected = painted.obs["louvain"].astype(str).map(palette).to_numpy()
    assert (drawn == expected).all()
    # the same colours come from a grammar plot with scale_color_obs
    grammar = ag.gganndata(painted, aes("UMAP_1", "UMAP_2", color="louvain")) + geom_point()
    grammar = grammar + ag.scale_color_obs(painted, "louvain")
    grammar._build()
    assert (_layer_data(grammar)["color"].to_numpy() == expected).all()


def test_filter_keeps_used_levels_in_order_and_colours_keyed_by_label(tidy):
    original = ag.obs_colors(tidy, "louvain")
    kept = tidy.filter(obs.louvain.is_in(["3", "1"]))
    assert list(kept.obs["louvain"].cat.categories) == ["1", "3"]  # unused dropped, order kept
    assert ag.obs_colors(kept, "louvain") == {"1": original["1"], "3": original["3"]}
    plot = ag.plot_embedding(kept, color="louvain")
    plot._build()
    expected = kept.obs["louvain"].astype(str).map(original).to_numpy()
    assert (_layer_data(plot)["color"].to_numpy() == expected).all()


# ---- count and summarize key order ------------------------------------------


def test_count_feeds_proportions_with_matching_numbers(tidy):
    counts = tidy.count(obs.phase, obs.louvain)
    assert list(counts.columns) == ["phase", "louvain", "n"]  # groups in key order
    plot = ag.plot_proportions(tidy, group_by="louvain", split_by="phase")
    merged = plot.data.merge(counts, on=["phase", "louvain"], suffixes=("", "_ap"), how="left")
    assert (merged["n"] == merged["n_ap"].fillna(0)).all()
    totals = plot.data.groupby("phase", observed=True)["value"].sum()
    np.testing.assert_allclose(totals, 1.0)


def test_summarize_feeds_a_bar_plot_in_level_order(tidy):
    summary = tidy.summarize(n=ap.n(), by=obs.louvain)
    assert summary["n"].sum() == tidy.n_obs
    plot = ggplot(summary, aes("louvain", "n")) + geom_col() + ag.theme_ggann()
    plot._build()
    levels = list(tidy.obs["louvain"].cat.categories)
    assert list(summary["louvain"].astype(str)) == levels
    np.testing.assert_array_equal(_layer_data(plot)["y"], summary["n"])


def test_input_is_not_mutated_by_plotting(adata, tidy):
    before = (adata.n_obs, adata.obs.columns.tolist(), adata.uns["louvain_colors"].copy())
    ag.plot_embedding(tidy.filter(obs.n_counts > 3333), color="louvain")._build()
    after = (adata.n_obs, adata.obs.columns.tolist(), adata.uns["louvain_colors"])
    assert before[:2] == after[:2]
    np.testing.assert_array_equal(before[2], after[2])
