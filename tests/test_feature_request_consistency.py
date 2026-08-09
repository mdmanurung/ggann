from __future__ import annotations

import importlib

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

import ggann as ag


def _feature_order(plot) -> list[str]:
    return list(plot.data["feature"].cat.categories)


def test_public_expression_helpers_deduplicate_features(adata, markers, group_key):
    genes = [markers[0], markers[0], markers[1]]
    expected = markers[:2]
    plots = [
        ag.plot_features(adata, genes, use_raw=False),
        ag.plot_dotplot(adata, genes, group_key, use_raw=False),
        ag.plot_matrixplot(adata, genes, group_key, use_raw=False),
        ag.plot_heatmap(adata, genes, group_key, use_raw=False),
        ag.plot_box(adata, genes, group_key, use_raw=False),
        ag.plot_sina(adata, genes, group_key, use_raw=False, bins=10),
        ag.plot_expression_bar(adata, genes, group_key, use_raw=False),
        ag.plot_expression_line(adata, genes, "phase", use_raw=False),
        ag.plot_violin(adata, genes, group_key, use_raw=False),
        ag.plot_stacked_violin(adata, genes, group_key, use_raw=False),
        ag.plot_ridge(adata, genes, group_key, use_raw=False, n_grid=16),
        ag.plot_tracksplot(adata, genes, group_key, use_raw=False),
    ]

    for plot in plots:
        order = _feature_order(plot)
        if plot is plots[3]:
            order = list(reversed(order))
        assert order == expected


def test_density_deduplicates_before_optional_backend(monkeypatch, adata, markers):
    module = importlib.import_module("ggann.density")
    calls: list[np.ndarray] = []

    def calculate_density(values, coordinates, **kwargs):
        calls.append(values)
        return np.arange(coordinates.shape[0], dtype=float)

    monkeypatch.setattr(module, "_require_pynebulosa", lambda: calculate_density)
    plot = ag.plot_density(
        adata,
        [markers[0], markers[0], markers[1]],
        use_raw=False,
    )

    assert len(calls) == 2
    assert _feature_order(plot) == markers[:2]


def test_correlation_counts_distinct_requested_genes(adata, markers, group_key):
    plot = ag.plot_correlation(
        adata,
        group_key,
        genes=[markers[0], markers[0], markers[1]],
        use_raw=False,
        cluster=False,
    )
    assert not plot.data.empty

    with pytest.raises(ValueError, match="at least two genes"):
        ag.plot_correlation(
            adata,
            group_key,
            genes=[markers[0], markers[0]],
            use_raw=False,
        )


def test_grouped_markers_deduplicate_within_one_group(adata, markers, group_key):
    plot = ag.plot_dotplot_grouped(
        adata,
        {"T cell": [markers[0], markers[0], markers[1]]},
        group_key,
        use_raw=False,
    )

    assert _feature_order(plot) == markers[:2]


def test_tidy_helper_accepts_nonunique_observation_names():
    with pytest.warns(UserWarning, match="Observation names are not unique"):
        adata = ad.AnnData(
            np.arange(6, dtype=float).reshape(3, 2),
            obs=pd.DataFrame(
                {"group": pd.Categorical(["a", "a", "b"])},
                index=["cell", "cell", "other"],
            ),
            var=pd.DataFrame(index=["g1", "g2"]),
        )

    plot = ag.plot_violin(adata, ["g1"], "group", use_raw=False)

    assert not isinstance(plot.data["obs_name"].dtype, pd.CategoricalDtype)
    assert set(plot.data["obs_name"]) == {"cell", "other"}


@pytest.mark.parametrize("groups", [["b", "a", "b"], ["a", "a", "b"]])
@pytest.mark.parametrize("helper", [ag.plot_heatmap, ag.plot_tracksplot])
def test_cell_rank_helpers_accept_nonunique_observation_names(helper, groups):
    with pytest.warns(UserWarning, match="Observation names are not unique"):
        adata = ad.AnnData(
            np.arange(6, dtype=float).reshape(3, 2),
            obs=pd.DataFrame(
                {"group": pd.Categorical(groups)},
                index=["cell", "cell", "other"],
            ),
            var=pd.DataFrame(index=["g1", "g2"]),
        )

    plot = helper(adata, ["g1", "g2"], "group", use_raw=False)

    assert plot.data["cell_rank"].notna().all()
    assert plot.data["cell_rank"].nunique() == adata.n_obs
    assert len(plot.data) == adata.n_obs * adata.n_vars


@pytest.mark.parametrize(
    "matrix_format",
    [
        pytest.param(lambda values: values, id="dense"),
        pytest.param(sparse.csr_matrix, id="csr"),
        pytest.param(sparse.csc_matrix, id="csc"),
    ],
)
def test_highest_expression_skips_nan_per_cell_and_gene(matrix_format):
    values = np.array(
        [
            [1.0, np.nan, 3.0],
            [2.0, 2.0, 0.0],
            [-0.5, 4.0, 1.0],
        ]
    )
    names = ["g1", "g2", "g3"]
    adata = ad.AnnData(
        matrix_format(values),
        obs=pd.DataFrame(index=["c1", "c2", "c3"]),
        var=pd.DataFrame(index=names),
    )
    wide = pd.DataFrame(values, columns=names)
    expected = wide.div(wide.sum(axis=1).replace(0, np.nan), axis=0) * 100.0
    expected_order = expected.mean().sort_values(ascending=False).index.tolist()

    with pytest.warns(UserWarning, match="negative values"):
        plot = ag.plot_highest_expr_genes(adata, n=3, use_raw=False)

    actual_order = list(reversed(plot.data["gene"].cat.categories))
    assert actual_order == expected_order
    for gene in names:
        actual = plot.data.loc[plot.data["gene"] == gene, "percent"].to_numpy()
        assert np.allclose(actual, expected[gene], equal_nan=True)


@pytest.mark.parametrize("source", ["x", "layer", "raw"])
@pytest.mark.parametrize(
    "matrix_format",
    [
        pytest.param(lambda values: values.copy(), id="dense"),
        pytest.param(sparse.csr_matrix, id="csr"),
        pytest.param(sparse.csc_matrix, id="csc"),
    ],
)
def test_highest_expression_zero_totals_match_scanpy_and_preserve_source(matrix_format, source):
    values = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 3.0, 0.0],
            [2.0, 2.0, 0.0],
            [9.0, 1.0, 0.0],
        ]
    )
    names = ["g1", "g2", "g3"]
    adata = ad.AnnData(
        matrix_format(values),
        obs=pd.DataFrame(index=["zero", "c2", "c3", "c4"]),
        var=pd.DataFrame(index=names),
    )
    adata.layers["counts"] = matrix_format(values)
    adata.raw = ad.AnnData(matrix_format(values), var=pd.DataFrame(index=names))
    kwargs = {
        "x": {"use_raw": False},
        "layer": {"layer": "counts"},
        "raw": {"use_raw": True},
    }[source]
    selected = {
        "x": adata.X,
        "layer": adata.layers["counts"],
        "raw": adata.raw.X,
    }[source]
    selected_before = selected.copy()

    with pytest.warns(UserWarning, match="zero total counts retained as 0%"):
        plot = ag.plot_highest_expr_genes(adata, n=3, **kwargs)

    denominators = values.sum(axis=1)
    denominators[denominators == 0] = 1.0
    expected = pd.DataFrame(values / denominators[:, None] * 100.0, columns=names)
    expected_order = expected.mean().sort_values(ascending=False, kind="stable").index.tolist()
    assert list(reversed(plot.data["gene"].cat.categories)) == expected_order
    for gene in names:
        actual = plot.data.loc[plot.data["gene"] == gene, "percent"].to_numpy()
        assert np.allclose(actual, expected[gene], rtol=1e-6, atol=1e-7, equal_nan=True)
        assert actual[0] == 0.0

    if sparse.issparse(selected_before):
        assert sparse.issparse(selected)
        assert selected.getformat() == selected_before.getformat()
        assert (selected != selected_before).nnz == 0
    else:
        np.testing.assert_array_equal(selected, selected_before)


@pytest.mark.parametrize(
    "matrix_format",
    [
        pytest.param(lambda values: values.copy(), id="dense"),
        pytest.param(sparse.csr_matrix, id="csr"),
        pytest.param(sparse.csc_matrix, id="csc"),
    ],
)
def test_highest_expression_zero_totals_support_backed_input(tmp_path, matrix_format):
    values = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 3.0, 0.0],
            [2.0, 2.0, 0.0],
        ],
        dtype=np.float32,
    )
    adata = ad.AnnData(
        matrix_format(values),
        var=pd.DataFrame(index=["g1", "g2", "g3"]),
    )
    path = tmp_path / "zero-totals.h5ad"
    adata.write_h5ad(path)
    backed = ad.read_h5ad(path, backed="r")
    matrix_type = type(backed.X)
    try:
        with pytest.warns(UserWarning, match="zero total counts retained as 0%"):
            plot = ag.plot_highest_expr_genes(backed, n=3, use_raw=False)
        assert backed.isbacked
        assert type(backed.X) is matrix_type
    finally:
        backed.file.close()

    for gene in ["g1", "g2", "g3"]:
        values_for_gene = plot.data.loc[plot.data["gene"] == gene, "percent"].to_numpy()
        assert values_for_gene[0] == 0.0


def test_highest_expression_supports_backed_dense_matrix(tmp_path):
    values = np.arange(1, 41, dtype=np.float32).reshape(10, 4)
    adata = ad.AnnData(
        values,
        obs=pd.DataFrame(index=[f"cell_{index}" for index in range(10)]),
        var=pd.DataFrame(index=[f"gene_{index}" for index in range(4)]),
    )
    expected = ag.plot_highest_expr_genes(adata, n=3, use_raw=False).data
    path = tmp_path / "dense.h5ad"
    adata.write_h5ad(path)
    backed = ad.read_h5ad(path, backed="r")
    try:
        actual = ag.plot_highest_expr_genes(backed, n=3, use_raw=False).data
    finally:
        backed.file.close()

    assert list(actual["gene"].cat.categories) == list(expected["gene"].cat.categories)
    assert np.allclose(actual["percent"], expected["percent"], equal_nan=True)


@pytest.mark.parametrize("n", [True, 1.5, "2"])
def test_highest_expression_rejects_non_integer_n(n):
    adata = ad.AnnData(
        np.ones((3, 2)),
        var=pd.DataFrame(index=["g1", "g2"]),
    )

    with pytest.raises(ValueError, match="positive integer"):
        ag.plot_highest_expr_genes(adata, n=n, use_raw=False)


def test_highest_expression_rejects_nonunique_variable_names():
    with pytest.warns(UserWarning, match="Variable names are not unique"):
        adata = ad.AnnData(
            np.arange(12, dtype=float).reshape(4, 3),
            var=pd.DataFrame(index=["g", "g", "h"]),
        )

    with pytest.raises(ValueError, match="Variable names.*unique"):
        ag.plot_highest_expr_genes(adata, n=3, use_raw=False)
