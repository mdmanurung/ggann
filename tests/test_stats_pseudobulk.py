"""Tests for the re-exported ggpubr stat layers and the pseudobulk wrapper."""

from __future__ import annotations

import numpy as np
import plotnine as p9
import pytest

import ggann as ag


def _build(plot):
    assert isinstance(plot, p9.ggplot)
    plot._build()
    return plot


STAT_LAYERS = [
    "stat_compare_means",
    "stat_pwc",
    "stat_pvalue_manual",
    "stat_cor",
    "stat_regline_equation",
    "stat_anova_test",
    "stat_kruskal_test",
    "stat_central_tendency",
    "geom_signif",
]


def test_stat_layers_reexported():
    for name in STAT_LAYERS:
        assert callable(getattr(ag, name)), name


def test_stat_cor_composes_on_scatter(adata):
    _build(ag.plot_qc_scatter(adata, x="n_counts", y="n_genes") + ag.stat_cor())


def test_stat_compare_means_composes_on_violin(adata, markers, group_key):
    _build(ag.plot_violin(adata, markers[:1], group_key) + ag.stat_central_tendency())


@pytest.fixture
def counts_adata(adata):
    ad = adata.copy()
    rng = np.random.RandomState(0)
    ad.obs["donor"] = rng.choice(["d1", "d2", "d3", "d4"], size=ad.n_obs)
    X = ad.raw[:, ad.var_names].X
    ad.layers["counts"] = np.rint(np.abs(X.toarray() if hasattr(X, "toarray") else np.asarray(X)))
    return ad


def test_pseudobulk_returns_anndata(counts_adata):
    pb = ag.pseudobulk(
        counts_adata,
        sample_col="donor",
        group_by="bulk_labels",
        layer="counts",
        min_cells=1,
        skip_checks=True,
    )
    from anndata import AnnData

    assert isinstance(pb, AnnData)
    assert pb.n_obs > 0
    assert {"donor", "bulk_labels"} <= set(pb.obs.columns)


def test_grammar_and_helpers_work_on_pseudobulk(counts_adata):
    # the whole point: the pseudobulk AnnData plugs straight back into ggann
    pb = ag.pseudobulk(
        counts_adata,
        sample_col="donor",
        group_by="bulk_labels",
        layer="counts",
        min_cells=1,
        skip_checks=True,
    )
    from plotnine import geom_boxplot

    from ggann import aes, gene, gganndata

    _build(gganndata(pb, aes("bulk_labels", gene("CD3D"), fill="donor")) + geom_boxplot())
    _build(ag.plot_dotplot(pb, ["CD3D", "NKG7", "CST3"], "bulk_labels", use_raw=False))
