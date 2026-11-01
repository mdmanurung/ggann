"""Guard the ggann side of the scanpy reproductions."""

from __future__ import annotations

import importlib.util
import os

import plotnine as p9
import pytest

_PATH = os.path.join(os.path.dirname(__file__), "..", "examples", "reproduce_scanpy.py")
_spec = importlib.util.spec_from_file_location("reproduce_scanpy", _PATH)
rs = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(rs)


@pytest.fixture(scope="module")
def de_adata(adata):
    import scanpy as sc

    ad = adata.copy()  # isolate: don't mutate the shared session fixture
    sc.tl.rank_genes_groups(ad, "bulk_labels", method="wilcoxon", n_genes=50, key_added=rs.DE_KEY)
    sc.tl.dendrogram(ad, groupby="bulk_labels")
    return ad


def test_ggann_reproductions_build(de_adata):
    figs = rs.pairs(de_adata)
    assert len(figs) >= 10
    for name, (_scfn, ggfn) in figs.items():
        plot = ggfn()
        assert isinstance(plot, p9.ggplot), name
        plot._build()
