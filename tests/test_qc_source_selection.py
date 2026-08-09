"""Regression tests for QC expression-source selection."""

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

import ggann as ag


def _qc_adata() -> AnnData:
    obs = pd.DataFrame(
        {"total_counts": [10, 20], "n_genes": [2, 3]},
        index=["cell_1", "cell_2"],
    )
    var = pd.DataFrame(index=["gene"])
    adata = AnnData(np.array([[1.0], [2.0]]), obs=obs, var=var)
    adata.layers["alternate"] = np.array([[11.0], [12.0]])
    adata.raw = AnnData(
        np.array([[21.0], [22.0]]),
        obs=obs.copy(),
        var=var.copy(),
    )
    return adata


def test_qc_scatter_selects_layer_for_gene_colour():
    adata = _qc_adata()

    plot = ag.plot_qc_scatter(
        adata,
        "total_counts",
        "n_genes",
        color="gene",
        layer="alternate",
    )

    np.testing.assert_array_equal(plot.data["gene"], [11.0, 12.0])


def test_qc_scatter_selects_raw_for_gene_colour():
    adata = _qc_adata()

    plot = ag.plot_qc_scatter(
        adata,
        "total_counts",
        "n_genes",
        color="gene",
        use_raw=True,
    )

    np.testing.assert_array_equal(plot.data["gene"], [21.0, 22.0])


@pytest.mark.parametrize(
    "x, y, color",
    [
        ("missing", "n_genes", None),
        ("total_counts", "missing", None),
        ("total_counts", "n_genes", "missing"),
    ],
)
def test_qc_scatter_rejects_unresolved_fields(x, y, color):
    with pytest.raises(KeyError, match="Could not resolve"):
        ag.plot_qc_scatter(_qc_adata(), x, y, color=color)
