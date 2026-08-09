"""Regression tests for the plot-aesthetics review fixes."""

from __future__ import annotations

import plotnine as p9
import pytest

import ggann as ag


def _build(plot):
    assert isinstance(plot, p9.ggplot)
    plot._build()
    return plot


def test_highest_expr_warns_on_scaled_matrix(adata):
    # pbmc68k_reduced.X is z-scored (has negatives) -> percentages are meaningless
    with pytest.warns(UserWarning, match="scaled"):
        ag.plot_highest_expr_genes(adata, n=5)


def test_highest_expr_on_raw_is_clean(adata):
    import warnings

    with warnings.catch_warnings():
        warnings.simplefilter("error")  # .raw (log-norm) must NOT trigger the scaled warning
        _build(ag.plot_highest_expr_genes(adata, n=10, use_raw=True))


def test_proportions_area_fills_missing_combos(adata, group_key):
    # drop one group entirely from one phase so a (phase, group) combo is missing;
    # the stacked area must still cover every split level (grid completed with 0).
    ad = adata.copy()
    victim = ad.obs[group_key].cat.categories[-1]
    keep = ~((ad.obs["phase"] == "G1") & (ad.obs[group_key] == victim)).to_numpy()
    ad = ad[keep].copy()
    _build(ag.plot_proportions(ad, group_key, split_by="phase", kind="area"))
    _build(ag.plot_proportions(ad, group_key, split_by="phase", kind="bar"))
