"""Smoke tests for the density / distribution / correlation helpers."""

from __future__ import annotations

import plotnine as p9
import pytest

import ggann as ag


def _is_ggplot(obj):
    assert isinstance(obj, p9.ggplot)
    # force the build so aesthetic / facet errors surface here, not only on save
    obj._build()
    return True


def test_plot_density_single(adata, markers):
    _is_ggplot(ag.plot_density(adata, markers[0], basis="umap"))


def test_plot_density_multi_joint(adata, markers):
    plot = ag.plot_density(adata, markers[:2], basis="umap", joint=True)
    _is_ggplot(plot)
    # one panel per feature plus the joint panel
    labels = list(plot.data["feature"].cat.categories)
    assert labels[-1] == " + ".join(markers[:2])
    assert len(labels) == 3


def test_plot_density_obs_feature(adata):
    _is_ggplot(ag.plot_density(adata, "n_genes", basis="umap"))


def test_plot_density_bad_feature_raises(adata):
    with pytest.raises(KeyError):
        ag.plot_density(adata, "not_a_real_gene_xyz", basis="umap")


def test_plot_box(adata, markers, group_key):
    _is_ggplot(ag.plot_box(adata, markers, group_key))


def test_plot_box_no_jitter(adata, markers, group_key):
    _is_ggplot(ag.plot_box(adata, markers, group_key, jitter=False))


def test_plot_expression_bar(adata, markers, group_key):
    _is_ggplot(ag.plot_expression_bar(adata, markers, group_key))


def test_plot_expression_bar_sd(adata, markers, group_key):
    _is_ggplot(ag.plot_expression_bar(adata, markers, group_key, error="sd"))


def test_plot_expression_line_categorical_x(adata, markers):
    _is_ggplot(ag.plot_expression_line(adata, markers, x="phase"))


def test_plot_expression_line_grouped(adata, markers, group_key):
    _is_ggplot(ag.plot_expression_line(adata, markers, x="phase", group_by=group_key))


def test_plot_correlation(adata, group_key):
    _is_ggplot(ag.plot_correlation(adata, group_key, cluster=True))


def test_plot_correlation_annotated_genes(adata, markers, group_key):
    _is_ggplot(ag.plot_correlation(adata, group_key, genes=markers, cluster=False, annotate=True))


def test_plot_embedding_labelled(adata, group_key):
    _is_ggplot(ag.plot_embedding(adata, "umap", color=group_key, label=True))


def test_repel_geoms_exported():
    # ggrepel-style non-overlapping labels, re-exported from plotnine-extra
    assert callable(ag.geom_text_repel)
    assert callable(ag.geom_label_repel)


def test_plot_upset(adata):
    genes = list(adata.raw.var_names[:30])
    upset = ag.plot_upset({"A": genes[:20], "B": genes[10:30], "C": genes[5:15]}, min_cardinality=1)
    assert upset.figure is not None


def test_plot_upset_needs_two_sets(adata):
    with pytest.raises(ValueError):
        ag.plot_upset({"only": ["a", "b"]})
