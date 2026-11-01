from __future__ import annotations

import importlib
import warnings
from pathlib import Path
from types import SimpleNamespace

import pytest

import ggann as ag


@pytest.fixture(scope="module")
def de_adata(adata):
    import scanpy as sc

    result = adata.copy()
    sc.tl.rank_genes_groups(result, "bulk_labels", method="wilcoxon", n_genes=20)
    return result


def test_rank_genes_helpers_accept_canonical_group_by(de_adata):
    with warnings.catch_warnings():
        warnings.simplefilter("error", FutureWarning)
        plot = ag.plot_rank_genes_matrixplot(de_adata, n_genes=2, group_by="bulk_labels")
    assert "bulk_labels" in plot.data.columns


def test_rank_genes_helpers_accept_scanpy_groupby(de_adata):
    # A1: groupby is a silent alias of group_by, as on every grouped helper
    with warnings.catch_warnings():
        warnings.simplefilter("error", FutureWarning)
        plot = ag.plot_rank_genes_dotplot(de_adata, n_genes=2, groupby="bulk_labels")
    assert "bulk_labels" in plot.data.columns


def test_rank_genes_helpers_reject_conflicting_group_spellings(de_adata):
    with pytest.raises(TypeError, match="only 'group_by'"):
        ag.plot_rank_genes_dotplot(
            de_adata,
            n_genes=2,
            group_by="bulk_labels",
            groupby="louvain",
        )


def test_pseudobulk_use_raw_and_legacy_alias(monkeypatch, adata):
    module = importlib.import_module("ggann.pseudobulk")
    calls = []

    def _pseudobulk(input_adata, **kwargs):
        calls.append({**kwargs, "input": input_adata})
        result = input_adata[:1].copy()
        result.obs["psbulk_cells"] = 10
        return result

    monkeypatch.setattr(
        module,
        "_require_decoupler",
        lambda: SimpleNamespace(pp=SimpleNamespace(pseudobulk=_pseudobulk)),
    )

    with warnings.catch_warnings():
        warnings.simplefilter("error", FutureWarning)
        ag.pseudobulk(adata, "phase", use_raw=True, min_cells=0)
    # raw is resolved to an AnnData holding the raw matrix, so decoupler gets raw=False
    assert calls[-1]["raw"] is False
    assert list(calls[-1]["input"].var_names) == list(adata.raw.var_names)
    with pytest.warns(FutureWarning, match="use_raw") as caught:
        ag.pseudobulk(adata, "phase", raw=True, min_cells=0)
    assert calls[-1]["raw"] is False
    assert Path(caught[0].filename) == Path(__file__)
    with pytest.raises(ValueError, match="layer"):
        ag.pseudobulk(adata, "phase", layer="counts", use_raw=True, min_cells=0)
