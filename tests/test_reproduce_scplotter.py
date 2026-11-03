"""Guard the scplotter reproduction figures + the plot_violin add_points option.

Building every reproduction here keeps the scplotter reproductions
honest -- if a helper or a composed layer stops producing a figure, this fails.
"""

from __future__ import annotations

import importlib.util
import os

import plotnine as p9
import pytest

import ggann as ag

_PATH = os.path.join(os.path.dirname(__file__), "..", "examples", "reproduce_scplotter.py")
_spec = importlib.util.spec_from_file_location("reproduce_scplotter", _PATH)
rs = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(rs)


@pytest.fixture(scope="module")
def de_adata(adata):
    import scanpy as sc

    if "rank_genes_groups" not in adata.uns:
        sc.tl.rank_genes_groups(adata, "bulk_labels", method="wilcoxon", n_genes=50)
    return adata


def test_scplotter_reproductions_build(de_adata):
    figs = rs.build(de_adata)
    assert len(figs) >= 15
    for name, (plot, _call) in figs.items():
        assert isinstance(plot, p9.ggplot), name
        plot._build()


def test_plot_violin_add_points(adata, markers, group_key):
    plot = ag.plot_violin(adata, markers[:2], group_key, add_points=True)
    assert isinstance(plot, p9.ggplot)
    plot._build()
