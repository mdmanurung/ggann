"""Values drawn by ggann helpers equal the values Scanpy's own helpers draw."""

import numpy as np
import pandas as pd
import pytest
import scanpy as sc

import ggann as ag

GENES = ["CD3D", "NKG7", "CST3"]


@pytest.fixture(scope="module")
def genes(adata):
    return [g for g in GENES if g in adata.raw.var_names]


def _scanpy_values(figure, adata):
    """Scanpy's plotted per-cell values as a long table keyed by cell and gene."""
    wide = figure.obs_tidy
    wide.index = adata.obs_names
    return wide.stack().rename("scanpy").rename_axis(["obs_name", "feature"])


def test_embedding_coordinates_equal_obsm(adata, group_key):
    plot = ag.plot_embedding(adata, "umap", color=group_key)
    assert np.allclose(plot.data[["UMAP_1", "UMAP_2"]].to_numpy(), adata.obsm["X_umap"][:, :2])


def test_violin_values_equal_scanpy_obs_tidy(adata, group_key, genes):
    fig = sc.pl.stacked_violin(
        adata, genes, groupby=group_key, use_raw=True, return_fig=True, show=False
    )
    ours = ag.plot_violin(adata, genes, group_key, use_raw=True).data
    ours = ours.set_index(["obs_name", "feature"])["value"]
    joined = pd.concat([ours, _scanpy_values(fig, adata)], axis=1, join="inner")
    assert len(joined) == adata.n_obs * len(genes)
    assert np.allclose(joined["value"], joined["scanpy"])


def test_stacked_violin_median_equals_group_median_of_scanpy_values(adata, group_key, genes):
    fig = sc.pl.stacked_violin(
        adata, genes, groupby=group_key, use_raw=True, return_fig=True, show=False
    )
    scanpy = _scanpy_values(fig, adata).reset_index()
    scanpy[group_key] = adata.obs[group_key].to_numpy()[
        adata.obs_names.get_indexer(scanpy["obs_name"])
    ]
    expected = scanpy.groupby([group_key, "feature"], observed=True)["scanpy"].median()
    ours = ag.plot_stacked_violin(adata, genes, group_key, use_raw=True).data
    ours = ours.groupby([group_key, "feature"], observed=True)["median"].first()
    joined = pd.concat([ours, expected], axis=1, join="inner")
    assert len(joined) == adata.obs[group_key].nunique() * len(genes)
    assert np.allclose(joined["median"], joined["scanpy"])


def test_heatmap_values_equal_scanpy_obs_tidy(adata, group_key, genes):
    fig = sc.pl.heatmap(adata, genes, groupby=group_key, use_raw=True, show=False, swap_axes=False)
    assert fig is not None  # Scanpy builds the same per-cell table the helper reads
    expected = sc.get.obs_df(adata, keys=genes, use_raw=True)
    ours = ag.plot_heatmap(adata, genes, group_key, use_raw=True).data
    wide = ours.pivot(index="obs_name", columns="feature", values="value").loc[
        expected.index, genes
    ]
    assert np.allclose(wide.to_numpy(), expected[genes].to_numpy())
