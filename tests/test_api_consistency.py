from __future__ import annotations

import importlib
from pathlib import Path
from types import SimpleNamespace
import warnings

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
        plot = ag.plot_rank_genes_matrixplot(
            de_adata, n_genes=2, group_by="bulk_labels"
        )
    assert "bulk_labels" in plot.data.columns


def test_rank_genes_helpers_warn_for_legacy_groupby(de_adata):
    with pytest.warns(FutureWarning, match="group_by") as caught:
        plot = ag.plot_rank_genes_dotplot(de_adata, n_genes=2, groupby="bulk_labels")
    assert "bulk_labels" in plot.data.columns
    assert Path(caught[0].filename) == Path(__file__)


def test_rank_genes_helpers_reject_both_group_spellings(de_adata):
    with pytest.raises(TypeError, match="only 'group_by'"):
        ag.plot_rank_genes_dotplot(
            de_adata,
            n_genes=2,
            group_by="bulk_labels",
            groupby="bulk_labels",
        )


def test_pseudobulk_use_raw_and_legacy_alias(monkeypatch, adata):
    module = importlib.import_module("ggann.pseudobulk")
    calls = []

    def _pseudobulk(input_adata, **kwargs):
        calls.append(kwargs)
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
    assert calls[-1]["raw"] is True
    with pytest.warns(FutureWarning, match="use_raw") as caught:
        ag.pseudobulk(adata, "phase", raw=True, min_cells=0)
    assert calls[-1]["raw"] is True
    assert Path(caught[0].filename) == Path(__file__)
    with pytest.raises(ValueError, match="layer"):
        ag.pseudobulk(adata, "phase", layer="counts", use_raw=True, min_cells=0)
