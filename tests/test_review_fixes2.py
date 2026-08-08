"""Regression tests for the code-review fixes."""

from __future__ import annotations

import numpy as np
import plotnine as p9
import pytest
import scipy.sparse as sp
from anndata import AnnData

import ggann as ag


def _build(plot):
    assert isinstance(plot, p9.ggplot)
    plot._build()
    return plot


def test_correlation_zero_variance_group_no_crash(adata, markers, group_key):
    # a group with a constant profile -> NaN correlations; clustering must not crash
    ad = adata.copy()
    X = ad.raw.X.toarray() if sp.issparse(ad.raw.X) else np.asarray(ad.raw.X)
    smallest = ad.obs[group_key].value_counts().index[-1]
    X[(ad.obs[group_key] == smallest).to_numpy(), :] = 5.0
    ad.raw = AnnData(X, var=ad.raw.var.copy())
    _build(
        ag.plot_correlation(
            ad, group_key, genes=list(ad.raw.var_names[:50]), cluster=True
        )
    )


def test_correlation_adaptive_and_forced_cmap(adata, markers, group_key):
    _build(
        ag.plot_correlation(adata, group_key, genes=markers)
    )  # cmap=None -> adaptive
    _build(ag.plot_correlation(adata, group_key, genes=markers, cmap="viridis"))


def test_correlation_two_groups_clusters(adata, markers):
    # exactly two groups must still honour cluster=True (guard is >= 2 now)
    ad = adata[adata.obs["bulk_labels"].isin(["CD19+ B", "Dendritic"])].copy()
    _build(ag.plot_correlation(ad, "bulk_labels", genes=markers, cluster=True))


def test_density_non_numeric_feature_raises(adata):
    with pytest.raises(TypeError, match="numeric"):
        ag.plot_density(adata, "phase", basis="umap")


def test_violin_box_and_points(adata, markers, group_key):
    _build(ag.plot_violin(adata, markers[:2], group_key, add_points=True))
    _build(
        ag.plot_violin(adata, markers[:1], group_key, add_box=False, add_points=True)
    )


def test_embedding_split_and_label(adata, group_key):
    _build(
        ag.plot_embedding(adata, "umap", color=group_key, split_by="phase", label=True)
    )


def test_embedding_label_warns_on_non_categorical(adata):
    with pytest.warns(UserWarning, match="label=True is ignored"):
        ag.plot_embedding(adata, "umap", color="CD3D", label=True)
    with pytest.warns(UserWarning, match="label=True is ignored"):
        ag.plot_embedding(adata, "umap", color=None, label=True)


def test_expression_line_name_collision_raises(adata, group_key):
    # a gene named the same as the group_by column would be dropped by melt -> guard
    with pytest.raises(ValueError, match="collide"):
        ag.plot_expression_line(adata, [group_key], x="phase", group_by=group_key)


def test_proportions_fill_labels_proportion(adata, group_key):
    plot = ag.plot_proportions(
        adata, group_key, split_by="phase", normalize=False, position="fill"
    )
    _build(plot)
    assert plot.labels.y == "proportion of cells"
