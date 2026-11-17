"""Scanpy keyword spellings, palette names and first-run error messages."""

from pathlib import Path

import pytest
import scanpy as sc

import ggann as ag

PBMC3K = Path(__file__).resolve().parents[1] / "data" / "pbmc3k_processed.h5ad"


def test_scanpy_spellings_give_the_same_data_as_ggann_spellings(adata, markers, group_key):
    ours = ag.plot_dotplot(adata, markers, group_key, cmap="Reds")
    theirs = ag.plot_dotplot(adata, var_names=markers, groupby=group_key, color_map="Reds")
    assert ours.data.equals(theirs.data)


def test_embedding_accepts_cmap_and_features_accept_color_map(adata, markers):
    ag.plot_embedding(adata, "umap", color=markers[0], cmap="magma")
    ag.plot_features(adata, markers[:2], color_map="magma")


def test_two_spellings_with_different_values_raise(adata, markers, group_key):
    with pytest.raises(TypeError, match="aliases"):
        ag.plot_dotplot(adata, markers, group_key, cmap="Reds", color_map="Blues")


def test_palette_accepts_a_colormap_name(adata, group_key):
    ag.plot_embedding(adata, "umap", color=group_key, palette="tab20")


def test_unknown_palette_name_raises_a_clear_error(adata, group_key):
    with pytest.raises(ValueError, match="not a Matplotlib colormap"):
        ag.plot_embedding(adata, "umap", color=group_key, palette="no_such_map")


def test_marker_dict_to_a_flat_list_helper_names_the_helpers(adata, markers, group_key):
    with pytest.raises(ValueError, match="plot_dotplot, plot_matrixplot"):
        ag.plot_violin(adata, {"T": markers[:2]}, group_key)


@pytest.mark.skipif(not PBMC3K.exists(), reason="shipped PBMC3K file not present")
def test_volcano_on_a_result_without_fold_changes_explains_the_fix():
    pbmc = sc.read_h5ad(PBMC3K)
    with pytest.raises(ValueError, match="re-run sc.tl.rank_genes_groups"):
        ag.plot_volcano(pbmc, group="B cells")


def test_volcano_unknown_group_lists_stored_groups(adata, group_key):
    local = adata.copy()
    sc.tl.rank_genes_groups(local, group_key, method="wilcoxon", key_added="w")
    with pytest.raises(ValueError, match="stored groups"):
        ag.plot_volcano(local, group="no such group", key="w")
    ag.plot_volcano(local, groups=local.obs[group_key].iloc[0], key="w")
