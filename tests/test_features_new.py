"""Tests for the feature-iteration additions: palette, DE, composition, faceting,
QC, marker summaries and helpers."""

import numpy as np
import pytest
from plotnine import ggplot

import ggann as ag
from ggann import obs_colors


def _renders(plot, tmp_path, name):
    out = tmp_path / f"{name}.png"
    plot.save(out, verbose=False, width=6, height=5, dpi=60)
    assert out.exists() and out.stat().st_size > 0


@pytest.fixture(scope="module")
def de_adata(adata):
    import scanpy as sc

    a = adata.copy()
    sc.tl.rank_genes_groups(a, "bulk_labels", method="wilcoxon", n_genes=50)
    return a


# --------------------------------------------------------------------------- #
# 1a. palette reuse
# --------------------------------------------------------------------------- #
def test_obs_colors_matches_uns(adata, group_key):
    mapping = obs_colors(adata, group_key)
    stored = [str(c) for c in adata.uns[f"{group_key}_colors"]]
    assert mapping is not None
    assert list(mapping.keys()) == list(adata.obs[group_key].cat.categories)
    assert list(mapping.values()) == stored


def test_obs_colors_none_for_noncategorical(adata):
    assert obs_colors(adata, "n_genes") is None


def test_embedding_categorical_renders_with_palette(adata, group_key, tmp_path):
    _renders(ag.plot_embedding(adata, "umap", color=group_key), tmp_path, "pal")


# --------------------------------------------------------------------------- #
# 1b. DE plots
# --------------------------------------------------------------------------- #
def test_rank_genes_df_columns(de_adata):
    df = ag.rank_genes_df(de_adata)
    assert {"group", "names", "logfoldchanges", "pvals_adj"} <= set(df.columns)


def test_de_requires_ranking(adata):
    a = adata.copy()
    a.uns.pop("rank_genes_groups", None)  # pbmc68k_reduced ships with one
    with pytest.raises(KeyError):
        ag.plot_rank_genes_dotplot(a)


def test_de_dotplot_and_volcano_render(de_adata, tmp_path):
    _renders(ag.plot_rank_genes_dotplot(de_adata, n_genes=3), tmp_path, "de_dot")
    _renders(ag.plot_rank_genes_matrixplot(de_adata, n_genes=3), tmp_path, "de_mat")
    _renders(ag.plot_volcano(de_adata, group="CD56+ NK"), tmp_path, "volcano")


# --------------------------------------------------------------------------- #
# 1c. composition
# --------------------------------------------------------------------------- #
def test_proportions_sum_to_one_per_split(adata, group_key):
    a = adata.copy()
    a.obs["cond"] = a.obs["phase"].astype(str)
    p = ag.plot_proportions(a, group_key, split_by="cond", normalize=True)
    sums = p.data.groupby("cond", observed=True)["value"].sum()
    assert np.allclose(sums.values, 1.0)


def test_proportions_counts_match_total(adata, group_key):
    p = ag.plot_proportions(adata, group_key, normalize=False)
    assert p.data["value"].sum() == adata.n_obs


# --------------------------------------------------------------------------- #
# 1d. faceted & split embeddings
# --------------------------------------------------------------------------- #
def test_split_embedding_renders(adata, group_key, tmp_path):
    a = adata.copy()
    a.obs["cond"] = a.obs["phase"].astype(str)
    _renders(ag.plot_embedding(a, "umap", color="CD3D", split_by="cond"), tmp_path, "split")


def test_plot_features_grid(adata, markers, tmp_path):
    p = ag.plot_features(adata, markers, basis="umap")
    assert isinstance(p, ggplot)
    assert set(p.data["feature"].unique()) == set(markers)
    _renders(p, tmp_path, "features")


def test_plot_features_needs_numeric(adata, group_key):
    with pytest.raises(ValueError, match="numeric feature"):
        ag.plot_features(adata, [group_key])  # categorical only


# --------------------------------------------------------------------------- #
# 2. QC + marker summaries
# --------------------------------------------------------------------------- #
def test_qc_violin_and_scatter(adata, group_key, tmp_path):
    _renders(ag.plot_qc_violin(adata, metrics=["n_genes", "n_counts"]), tmp_path, "qcv")
    _renders(ag.plot_qc_violin(adata, metrics=["n_genes"], group_by=group_key), tmp_path, "qcvg")
    _renders(ag.plot_qc_scatter(adata, "n_counts", "n_genes", color=group_key), tmp_path, "qcs")


def test_highest_expr_genes(adata, tmp_path):
    p = ag.plot_highest_expr_genes(adata, n=10)
    assert p.data["gene"].nunique() == 10
    _renders(p, tmp_path, "hexpr")


def test_stacked_violin_and_tracksplot(adata, markers, group_key, tmp_path):
    _renders(ag.plot_stacked_violin(adata, markers, group_key), tmp_path, "sv")
    _renders(ag.plot_tracksplot(adata, markers, group_key), tmp_path, "tp")


def test_gene_group_brackets(adata, group_key, tmp_path):
    gg = {"T": ["CD3D", "CD8A"], "NK": ["NKG7", "GNLY"]}
    p = ag.plot_dotplot_grouped(adata, gg, group_key)
    assert "gene_group" in p.data.columns
    _renders(p, tmp_path, "grouped")


# --------------------------------------------------------------------------- #
# 3. helpers
# --------------------------------------------------------------------------- #
def test_composition_reexports_present():
    for name in ["Beside", "Stack", "Wrap", "plot_layout", "plot_annotation", "ggsave"]:
        assert hasattr(ag, name)


def test_violin_stats_option(adata, markers, group_key, tmp_path):
    _renders(ag.plot_violin(adata, markers[:2], group_key, stats=True), tmp_path, "vstats")


# --------------------------------------------------------------------------- #
# regressions for review-confirmed defects
# --------------------------------------------------------------------------- #
def test_plot_features_dedup_duplicates(adata):
    p = ag.plot_features(adata, ["CD3D", "CD3D", "CD8A"])  # must not raise
    assert set(p.data["feature"].unique()) == {"CD3D", "CD8A"}


def test_split_by_equals_color(adata, group_key, tmp_path):
    # split_by and color naming the same obs column must not crash on the join
    _renders(
        ag.plot_embedding(adata, "umap", color=group_key, split_by=group_key),
        tmp_path,
        "splitcolor",
    )


def test_volcano_requires_pvalues(adata):
    # pbmc68k_reduced ships with a logreg rank_genes_groups (no pvals_adj)
    with pytest.raises(ValueError, match="missing"):
        ag.plot_volcano(adata, group="CD4+/CD25 T Reg")


def test_gene_in_two_groups_raises(adata, group_key):
    with pytest.raises(ValueError, match="more than one gene group"):
        ag.plot_dotplot_grouped(adata, {"A": ["CD3D"], "B": ["CD3D"]}, group_key)


def test_rank_genes_matrixplot(de_adata):
    _ = ag.plot_rank_genes_matrixplot(de_adata, n_genes=2)
