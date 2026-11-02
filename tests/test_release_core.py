"""Regressions for the 0.1.0 readiness rows N1/C5, B1, B2-B4 and B8.

Row IDs refer to plans/evaluation/remaining_todos.md.
"""

from __future__ import annotations

import warnings

import anndata as ad
import numpy as np
import pandas as pd
import pytest

import ggann as ag

GROUP = "bulk_labels"
GENES = ["CD3D", "NKG7"]


@pytest.fixture(scope="module")
def nan_adata(adata):
    """25 cells with a missing group and 10 more with a missing split."""
    a = adata.copy()
    a.obs.loc[a.obs_names[:25], GROUP] = np.nan
    a.obs.loc[a.obs_names[30:40], "depth"] = np.nan
    return a


GROUPED = {
    "box": lambda a, s: ag.plot_box(a, GENES, group_by=GROUP, split_by=s, use_raw=True),
    "violin": lambda a, s: ag.plot_violin(a, GENES, group_by=GROUP, split_by=s, use_raw=True),
    "sina": lambda a, s: ag.plot_sina(a, GENES, group_by=GROUP, split_by=s, use_raw=True),
    "exprbar": lambda a, s: ag.plot_expression_bar(
        a, GENES, group_by=GROUP, split_by=s, use_raw=True
    ),
    "dotplot": lambda a, s: ag.plot_dotplot(a, GENES, group_by=GROUP, split_by=s, use_raw=True),
    "matrix": lambda a, s: ag.plot_matrixplot(a, GENES, group_by=GROUP, split_by=s, use_raw=True),
    "proportions": lambda a, s: ag.plot_proportions(a, group_by=GROUP, split_by=s),
}
UNSPLIT = {
    "ridge": lambda a: ag.plot_ridge(a, GENES, group_by=GROUP, use_raw=True),
    "stacked": lambda a: ag.plot_stacked_violin(a, GENES, group_by=GROUP, use_raw=True),
    "heatmap": lambda a: ag.plot_heatmap(a, GENES, group_by=GROUP, use_raw=True),
    "tracksplot": lambda a: ag.plot_tracksplot(a, GENES, group_by=GROUP, use_raw=True),
    "qc_violin": lambda a: ag.plot_qc_violin(a, ["n_counts"], group_by=GROUP),
    "embedding_density": lambda a: ag.plot_embedding_density(a, "umap", group_by=GROUP),
}


@pytest.mark.parametrize("name", GROUPED)
def test_missing_group_and_split_dropped_with_count(nan_adata, name):
    with pytest.warns(
        UserWarning, match=r"dropped 35 of 700 cells with missing 'bulk_labels', 'depth'"
    ):
        plot = GROUPED[name](nan_adata, "depth")
    plot.draw()
    if GROUP in plot.data:
        assert not plot.data[GROUP].isna().any()


@pytest.mark.parametrize("name", UNSPLIT)
def test_missing_group_dropped_and_draws(nan_adata, name):
    with pytest.warns(UserWarning, match=r"dropped 25 of 700 cells with missing 'bulk_labels'"):
        plot = UNSPLIT[name](nan_adata)
    plot.draw()
    if GROUP in plot.data:
        assert not plot.data[GROUP].isna().any()


def test_no_missing_groups_no_warning(adata):
    with warnings.catch_warnings():
        warnings.filterwarnings("error", message="dropped .* cells")
        ag.plot_violin(adata, GENES, group_by=GROUP, use_raw=True)


def test_proportions_use_kept_cells_as_denominator(nan_adata):
    with pytest.warns(UserWarning, match="dropped"):
        plot = ag.plot_proportions(nan_adata, group_by=GROUP, split_by="depth", normalize=True)
    # each split level's composition sums to 1 over the kept cells only
    data = plot.data
    totals = data.groupby("depth", observed=True)["value"].sum()
    np.testing.assert_allclose(totals.to_numpy(), 1.0)
    kept = nan_adata.obs.dropna(subset=[GROUP, "depth"])
    assert int(data["n"].sum()) == len(kept)


def test_input_not_mutated_by_drop(nan_adata):
    before = nan_adata.obs[GROUP].isna().sum()
    with pytest.warns(UserWarning):
        ag.plot_violin(nan_adata, GENES, group_by=GROUP, use_raw=True)
    assert nan_adata.obs[GROUP].isna().sum() == before == 25


# B1: a bare string is one name; an empty request fails at the call.
@pytest.mark.parametrize(
    "call",
    [
        lambda a, g: ag.plot_violin(a, g, group_by=GROUP, use_raw=True),
        lambda a, g: ag.plot_dotplot(a, g, group_by=GROUP, use_raw=True),
        lambda a, g: ag.plot_ridge(a, g, group_by=GROUP, use_raw=True),
        lambda a, g: ag.plot_features(a, g, use_raw=True),
    ],
)
def test_scalar_gene_string_is_one_name(adata, call):
    plot = call(adata, "CD3D")
    assert list(map(str, plot.data["feature"].unique())) == ["CD3D"]


@pytest.mark.parametrize(
    "call",
    [
        lambda a: ag.plot_violin(a, [], group_by=GROUP),
        lambda a: ag.plot_ridge(a, [], group_by=GROUP),
        lambda a: ag.plot_dotplot(a, [], group_by=GROUP),
        lambda a: ag.plot_expression_line(a, [], x="n_counts"),
        lambda a: ag.plot_qc_violin(a, [], group_by=GROUP),
    ],
)
def test_empty_request_raises(adata, call):
    with pytest.raises(ValueError, match="No (genes|QC metrics) requested"):
        call(adata)


def test_qc_metric_string_is_one_name(adata):
    plot = ag.plot_qc_violin(adata, "n_genes", group_by=GROUP)
    plot.draw()


# B2-B4: helper fields are validated at the call, never silently skipped.
def test_missing_obs_column_raises_at_call(adata):
    with pytest.raises(KeyError, match="'zzz'"):
        ag.plot_violin(adata, ["CD3D"], group_by=GROUP, split_by="zzz", use_raw=True)
    with pytest.raises(KeyError, match="'zzz'"):
        ag.plot_embedding(adata, "umap", color=GROUP, split_by="zzz")


def test_features_reject_unknown_and_categorical(adata):
    with pytest.raises(KeyError, match="'NOPE'"):
        ag.plot_features(adata, ["CD3D", "NOPE"], use_raw=True)
    with pytest.raises(ValueError, match="categorical"):
        ag.plot_features(adata, ["CD3D", "louvain"], use_raw=True)


def test_bool_color_is_discrete(adata):
    a = adata.copy()
    a.obs["flag"] = a.obs["n_counts"] > a.obs["n_counts"].median()
    ag.plot_embedding(a, "umap", color="flag").draw()


# B8: the other expression source is named in the error.
def test_gene_in_other_source_hint():
    rng = np.random.default_rng(0)
    full = ad.AnnData(
        rng.poisson(1, (40, 6)).astype(float),
        obs=pd.DataFrame(
            {"g": pd.Categorical(rng.choice(["a", "b"], 40))},
            index=[f"c{i}" for i in range(40)],
        ),
        var=pd.DataFrame(index=[f"G{i}" for i in range(6)]),
    )
    full.raw = full
    a = full[:, :3].copy()
    with pytest.raises(KeyError, match=r"They are in \.raw; pass use_raw=True"):
        ag.plot_violin(a, ["G5"], group_by="g", use_raw=False)


# A1: Scanpy's groupby= spelling is a silent alias on every group_by helper.
ALIAS_CALLS = {
    "violin": lambda a, **k: ag.plot_violin(a, GENES, use_raw=True, **k),
    "dotplot": lambda a, **k: ag.plot_dotplot(a, GENES, use_raw=True, **k),
    "matrixplot": lambda a, **k: ag.plot_matrixplot(a, GENES, use_raw=True, **k),
    "heatmap": lambda a, **k: ag.plot_heatmap(a, GENES, use_raw=True, **k),
    "tracksplot": lambda a, **k: ag.plot_tracksplot(a, GENES, use_raw=True, **k),
    "ridge": lambda a, **k: ag.plot_ridge(a, GENES, use_raw=True, **k),
    "proportions": lambda a, **k: ag.plot_proportions(a, **k),
    "qc_violin": lambda a, **k: ag.plot_qc_violin(a, ["n_counts"], **k),
}


@pytest.mark.parametrize("name", ALIAS_CALLS)
def test_A1_groupby_alias_matches_group_by(adata, name):
    call = ALIAS_CALLS[name]
    expected = call(adata, group_by=GROUP).data
    got = call(adata, groupby=GROUP).data
    pd.testing.assert_frame_equal(got.reset_index(drop=True), expected.reset_index(drop=True))
    with pytest.raises(TypeError, match="aliases"):
        call(adata, group_by=GROUP, groupby="louvain")


def test_A1_every_group_by_helper_is_aliased():
    import inspect

    missing = [
        name
        for name in ag.__all__
        if inspect.isfunction(getattr(ag, name))
        and "group_by" in inspect.signature(getattr(ag, name)).parameters
        and not hasattr(getattr(ag, name), "__wrapped__")
    ]
    assert missing == []


def test_warnings_point_at_the_callers_line(nan_adata):
    # skip_file_prefixes: decorated helpers still attribute warnings to user code
    from pathlib import Path

    with pytest.warns(UserWarning, match="dropped") as caught:
        ag.plot_violin(nan_adata, GENES, groupby=GROUP, use_raw=True)
    assert Path(caught[0].filename) == Path(__file__)
