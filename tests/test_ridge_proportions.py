"""Tests for plot_ridge and the plot_proportions kind= variants."""

from __future__ import annotations

import numpy as np
import pandas as pd
import plotnine as p9
import pytest
import scipy.sparse as sp
from anndata import AnnData

import ggann as ag


def _build(plot):
    assert isinstance(plot, p9.ggplot)
    plot._build()
    return plot


def test_plot_ridge(adata, markers, group_key):
    _build(ag.plot_ridge(adata, markers[:2], group_key))
    _build(ag.plot_ridge(adata, markers[:1], group_key))


def test_plot_ridge_degenerate_group(adata, markers, group_key):
    # a constant / tiny group must draw a flat baseline, not crash the KDE
    ad = adata.copy()
    X = ad.raw.X.toarray() if sp.issparse(ad.raw.X) else np.asarray(ad.raw.X)
    smallest = ad.obs[group_key].value_counts().index[-1]
    X[(ad.obs[group_key] == smallest).to_numpy(), :] = 0.0
    ad.raw = AnnData(X, var=ad.raw.var.copy())
    _build(ag.plot_ridge(ad, markers[:1], group_key))


def test_proportions_kinds(adata, group_key):
    _build(ag.plot_proportions(adata, group_key, split_by="phase", kind="area"))
    _build(ag.plot_proportions(adata, group_key, split_by="phase", kind="trend"))
    _build(ag.plot_proportions(adata, group_key, split_by="phase", kind="bar"))


def test_proportions_area_trend_need_split(adata, group_key):
    with pytest.raises(ValueError, match="split_by"):
        ag.plot_proportions(adata, group_key, kind="area")
    with pytest.raises(ValueError, match="split_by"):
        ag.plot_proportions(adata, group_key, kind="trend")


def test_proportions_bad_kind(adata, group_key):
    with pytest.raises(ValueError, match="kind"):
        ag.plot_proportions(adata, group_key, kind="pie")


@pytest.mark.parametrize(
    "group_by, split_by", [("missing", None), ("group", "missing")]
)
def test_proportions_rejects_missing_grouping_columns(group_by, split_by):
    adata = AnnData(
        np.ones((2, 1)),
        obs=pd.DataFrame({"group": ["a", "b"]}, index=["cell_1", "cell_2"]),
        var=pd.DataFrame(index=["gene"]),
    )

    with pytest.raises(KeyError, match="Observation column"):
        ag.plot_proportions(adata, group_by, split_by=split_by)


def test_proportions_drops_unused_categorical_levels():
    adata = AnnData(
        np.ones((4, 1)),
        obs=pd.DataFrame(
            {
                "group": pd.Categorical(
                    ["a", "a", "b", "b"], categories=["a", "b", "empty"]
                ),
                "split": pd.Categorical(
                    ["x", "y", "x", "y"], categories=["x", "y", "empty"]
                ),
            },
            index=[f"cell_{index}" for index in range(4)],
        ),
        var=pd.DataFrame(index=["gene"]),
    )

    plot = ag.plot_proportions(adata, "group", split_by="split")

    assert list(plot.data["group"].cat.categories) == ["a", "b"]
    assert list(plot.data["split"].cat.categories) == ["x", "y"]
    assert plot.data["value"].notna().all()
