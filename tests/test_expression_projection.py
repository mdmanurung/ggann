from __future__ import annotations

import anndata as ad
import annplyr as ap
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

import ggann as ag
from ggann._aggregate import aggregate_expression, group_means, tidy_expression
from ggann._expression import (
    densify_frame,
    expression_frame,
    project_expression,
    resolve_source,
)
from ggann._resolve import resolve_frame


def _adata(matrix_format: str = "dense") -> ad.AnnData:
    rng = np.random.default_rng(4)
    values = rng.poisson(0.4, size=(40, 24)).astype(np.float32)
    if matrix_format == "csr":
        matrix = sparse.csr_matrix(values)
    elif matrix_format == "csc":
        matrix = sparse.csc_matrix(values)
    else:
        matrix = values
    obs = pd.DataFrame(
        {
            "group": pd.Categorical(np.tile(["a", "b", "c", "d"], 10)),
            "batch": pd.Categorical(np.repeat(["x", "y"], 20)),
        },
        index=[f"cell_{index}" for index in range(values.shape[0])],
    )
    result = ad.AnnData(
        X=matrix,
        obs=obs,
        var=pd.DataFrame(index=[f"gene_{index}" for index in range(values.shape[1])]),
    )
    result.layers["double"] = matrix * 2
    result.raw = result.copy()
    result.obsm["X_pca"] = rng.normal(size=(result.n_obs, 50)).astype(np.float32)
    result.obsp["graph"] = sparse.eye(result.n_obs, format="csr")
    result.uns["metadata"] = {"large": np.ones(1_000)}
    return result


@pytest.mark.parametrize("matrix_format", ["dense", "csr", "csc"])
def test_projection_limits_matrix_before_annplyr(matrix_format):
    adata = _adata(matrix_format)
    projected, genes = project_expression(
        adata,
        ["gene_7", "gene_2", "gene_7"],
        kind="x",
        layer=None,
        obs=["group"],
    )

    assert genes == ["gene_7", "gene_2"]
    assert projected.shape == (adata.n_obs, 2)
    assert list(projected.obs.columns) == ["group"]
    assert not projected.obsm
    assert not projected.obsp
    assert not projected.uns
    expected = np.asarray(
        adata[:, genes].X.toarray() if sparse.issparse(adata.X) else adata[:, genes].X
    )
    actual = projected.X.toarray() if sparse.issparse(projected.X) else np.asarray(projected.X)
    assert np.array_equal(actual, expected)


def test_raw_projection_does_not_copy_graphs_or_uns():
    adata = _adata("csr")
    projected, genes = project_expression(
        adata,
        ["gene_3", "gene_9"],
        kind="raw",
        layer=None,
        obs=["batch"],
    )

    assert genes == ["gene_3", "gene_9"]
    assert projected.shape == (adata.n_obs, 2)
    assert not projected.obsm
    assert not projected.obsp
    assert not projected.uns
    assert "graph" in adata.obsp
    assert "metadata" in adata.uns


@pytest.mark.parametrize("matrix_format", ["dense", "csr", "csc"])
@pytest.mark.parametrize(
    ("kind", "layer", "multiplier"),
    [("x", None, 1.0), ("layer", "double", 2.0), ("raw", None, 1.0)],
)
def test_expression_frame_uses_projected_public_sources(matrix_format, kind, layer, multiplier):
    adata = _adata(matrix_format)
    result, genes = expression_frame(
        adata,
        ["gene_9", "gene_2", "gene_9"],
        kind=kind,
        layer=layer,
        max_matrix_values=adata.n_obs * 2,
    )

    assert genes == ["gene_9", "gene_2"]
    assert list(result.columns) == genes
    assert result.index.equals(adata.obs_names)
    matrix = adata.X.toarray() if sparse.issparse(adata.X) else np.asarray(adata.X)
    expected = matrix[:, [9, 2]] * multiplier
    actual = densify_frame(result).to_numpy()
    assert np.array_equal(actual, expected)


def test_expression_frame_propagates_annplyr_budget_error():
    adata = _adata("csr")
    with pytest.raises(ap.AnnplyrError, match="materialize 80 matrix values"):
        expression_frame(
            adata,
            ["gene_1", "gene_2"],
            kind="x",
            layer=None,
            max_matrix_values=79,
        )


@pytest.mark.parametrize("matrix_format", ["dense", "csr", "csc"])
def test_backed_raw_resolution_and_aggregation_match_memory(matrix_format, tmp_path):
    adata = _adata(matrix_format)
    fields = [ag.gene("gene_3", use_raw=True), "group"]
    expected_frame = resolve_frame(adata, fields, use_raw=False)
    expected_aggregate = aggregate_expression(adata, ["gene_3", "gene_9"], "group", use_raw=True)
    path = tmp_path / f"raw-{matrix_format}.h5ad"
    adata.write_h5ad(path)
    backed = ad.read_h5ad(path, backed="r")
    try:
        actual_frame = resolve_frame(backed, fields, use_raw=False)
        actual_aggregate = aggregate_expression(backed, ["gene_3", "gene_9"], "group", use_raw=True)
    finally:
        backed.file.close()

    pd.testing.assert_frame_equal(expected_frame, actual_frame, check_dtype=False)
    pd.testing.assert_frame_equal(
        expected_aggregate,
        actual_aggregate,
        check_dtype=False,
        check_exact=False,
        rtol=1e-6,
        atol=1e-7,
    )


@pytest.mark.parametrize("matrix_format", ["dense", "csr", "csc"])
def test_backed_mixed_sources_and_reordered_view_match_memory(matrix_format, tmp_path):
    adata = _adata(matrix_format)
    positions = [11, 2, 29, 0, 17]
    fields = [
        ag.gene("gene_3", use_raw=False),
        ag.gene("gene_9", layer="double"),
        ag.gene("gene_7", use_raw=True),
        ag.obsm("pca", 12),
        "group",
    ]
    expected = resolve_frame(adata[positions], fields, use_raw=False)
    path = tmp_path / f"mixed-{matrix_format}.h5ad"
    adata.write_h5ad(path)
    backed = ad.read_h5ad(path, backed="r")
    try:
        actual = resolve_frame(backed[positions], fields, use_raw=False)
    finally:
        backed.file.close()

    pd.testing.assert_frame_equal(expected, actual, check_dtype=False)


@pytest.mark.parametrize("matrix_format", ["dense", "csr", "csc"])
@pytest.mark.parametrize("expression_cutoff", [0.0, -0.5])
def test_projected_aggregation_matches_manual_means_and_fractions(matrix_format, expression_cutoff):
    adata = _adata(matrix_format)
    genes = ["gene_5", "gene_1", "gene_5", "gene_8"]
    result = aggregate_expression(
        adata,
        genes,
        "group",
        use_raw=False,
        expression_cutoff=expression_cutoff,
    )

    matrix = adata.X.toarray() if sparse.issparse(adata.X) else np.asarray(adata.X)
    wide = pd.DataFrame(matrix, index=adata.obs_names, columns=adata.var_names)
    wide["group"] = adata.obs["group"]
    expected_genes = ["gene_5", "gene_1", "gene_8"]
    expected_mean = wide.groupby("group", observed=True)[expected_genes].mean()
    expected_fraction = (
        wide.assign(**{gene: wide[gene] > expression_cutoff for gene in expected_genes})
        .groupby("group", observed=True)[expected_genes]
        .mean()
    )

    actual_mean = result.pivot(index="group", columns="feature", values="mean_expression")
    actual_fraction = result.pivot(index="group", columns="feature", values="fraction")
    assert list(result["feature"].cat.categories) == expected_genes
    assert np.allclose(actual_mean.loc[:, expected_genes], expected_mean)
    assert np.allclose(actual_fraction.loc[:, expected_genes], expected_fraction)


@pytest.mark.parametrize("matrix_format", ["dense", "csr", "csc"])
@pytest.mark.parametrize("expression_cutoff", [0.0, -0.5])
def test_projected_aggregation_handles_expression_nan_like_pandas(matrix_format, expression_cutoff):
    values = np.array(
        [[np.nan, 0.0], [1.0, -1.0], [0.0, 2.0], [3.0, np.nan]],
        dtype=np.float32,
    )
    matrix = (
        sparse.csr_matrix(values)
        if matrix_format == "csr"
        else sparse.csc_matrix(values)
        if matrix_format == "csc"
        else values
    )
    adata = ad.AnnData(
        X=matrix,
        obs=pd.DataFrame(
            {"group": pd.Categorical(["a", "a", "b", "b"])},
            index=[f"cell_{index}" for index in range(4)],
        ),
        var=pd.DataFrame(index=["gene_0", "gene_1"]),
    )

    result = aggregate_expression(
        adata,
        ["gene_0", "gene_1"],
        "group",
        use_raw=False,
        expression_cutoff=expression_cutoff,
    )
    wide = pd.DataFrame(values, columns=adata.var_names)
    wide["group"] = adata.obs["group"].to_numpy()
    expected_mean = wide.groupby("group", observed=True)[list(adata.var_names)].mean()
    expected_fraction = (
        wide.assign(
            **{gene: wide[gene].gt(expression_cutoff).fillna(False) for gene in adata.var_names}
        )
        .groupby("group", observed=True)[list(adata.var_names)]
        .mean()
    )
    actual_mean = result.pivot(index="group", columns="feature", values="mean_expression")
    actual_fraction = result.pivot(index="group", columns="feature", values="fraction")
    assert np.allclose(actual_mean, expected_mean, equal_nan=True)
    assert np.allclose(actual_fraction, expected_fraction, equal_nan=True)


def test_tidy_expression_deduplicates_requested_genes():
    adata = _adata("csr")
    result = tidy_expression(adata, ["gene_1", "gene_1", "gene_2"], "group", use_raw=False)
    assert list(result["feature"].cat.categories) == ["gene_1", "gene_2"]
    assert len(result) == adata.n_obs * 2


def test_aggregation_accepts_nonunique_observation_names():
    with pytest.warns(UserWarning, match="Observation names are not unique"):
        adata = ad.AnnData(
            np.array([[1.0], [10.0], [100.0]]),
            obs=pd.DataFrame(
                {"group": pd.Categorical(["b", "a", "b"])},
                index=["cell", "cell", "other"],
            ),
            var=pd.DataFrame(index=["gene"]),
        )

    result = aggregate_expression(adata, ["gene"], "group", use_raw=False)
    mean = result.set_index("group")["mean_expression"]
    fraction = result.set_index("group")["fraction"]

    assert len(result) == 2
    assert mean["a"] == 10.0
    assert mean["b"] == 50.5
    assert fraction.to_dict() == {"a": 1.0, "b": 1.0}


def test_mixed_resolution_is_positional_with_duplicate_reordered_obs_names():
    adata = _adata("csr")
    adata.obs_names = ["cell" if index % 3 else "other" for index in range(40)]
    positions = [9, 2, 8, 1, 7]
    view = adata[positions]
    fields = [
        ag.gene("gene_4", use_raw=False),
        ag.gene("gene_6", layer="double"),
        ag.obsm("pca", 3),
        "batch",
    ]

    result = resolve_frame(view, fields, use_raw=False)

    assert result.index.tolist() == view.obs_names.tolist()
    assert result["batch"].tolist() == view.obs["batch"].tolist()
    assert np.array_equal(result["gene_4"].to_numpy(), adata.X[positions, 4].toarray().ravel())
    assert np.array_equal(
        result["gene_6"].to_numpy(),
        adata.layers["double"][positions, 6].toarray().ravel(),
    )
    assert np.array_equal(result["PC_4"].to_numpy(), adata.obsm["X_pca"][positions, 3])


def test_resolution_preserves_categorical_order_and_missing_values():
    adata = _adata("csr")
    values = np.tile(["a", "b", None, "c"], 10)
    adata.obs["group"] = pd.Categorical(
        values,
        categories=["c", "b", "a", "unused"],
        ordered=True,
    )

    result = resolve_frame(
        adata,
        ["group", ag.gene("gene_2", use_raw=False)],
        use_raw=False,
    )

    assert isinstance(result["group"].dtype, pd.CategoricalDtype)
    assert result["group"].cat.categories.tolist() == ["c", "b", "a", "unused"]
    assert result["group"].cat.ordered
    assert result["group"].isna().sum() == 10


def test_resolve_frame_preflights_cumulative_budget_before_any_accessor_call(
    monkeypatch,
):
    adata = _adata("csr")
    accessor_type = type(adata.ap)
    original = accessor_type.to_df
    calls = []

    def recording_to_df(self, *args, **kwargs):
        calls.append(kwargs.copy())
        return original(self, *args, **kwargs)

    monkeypatch.setattr(accessor_type, "to_df", recording_to_df)
    fields = [
        ag.gene("gene_1", use_raw=False),
        ag.gene("gene_2", layer="double"),
        ag.obsm("pca", 5),
        "group",
    ]
    projected_values = adata.n_obs * 3

    with pytest.raises(ap.AnnplyrError, match=f"materialize {projected_values}"):
        resolve_frame(
            adata,
            fields,
            use_raw=False,
            max_matrix_values=projected_values - 1,
        )
    assert calls == []

    result = resolve_frame(
        adata,
        fields,
        use_raw=False,
        max_matrix_values=projected_values,
    )
    assert list(result.columns) == ["gene_1", "gene_2", "PC_6", "group"]
    assert len(calls) == 2
    assert calls[0]["obs"] == ["group"]
    assert calls[0]["obsm"] == {"X_pca": ["5"]}
    assert calls[0]["x"] == ["gene_1"]
    assert calls[1]["layer"] == "double"


def test_grouped_aggregation_fuses_expression_and_obs_projection(monkeypatch):
    adata = _adata("csr")
    accessor_type = type(adata.ap)
    original = accessor_type.to_df
    calls = []

    def recording_to_df(self, *args, **kwargs):
        calls.append(kwargs.copy())
        return original(self, *args, **kwargs)

    monkeypatch.setattr(accessor_type, "to_df", recording_to_df)
    result = aggregate_expression(adata, ["gene_1", "gene_3"], "group", use_raw=False)

    assert not result.empty
    assert len(calls) == 1
    assert calls[0]["obs"] == ["group"]
    assert calls[0]["x"] == ["gene_1", "gene_3"]


def test_gganndata_exposes_cumulative_materialization_budget():
    adata = _adata("csr")
    mapping = ag.aes(
        x=ag.obsm("pca", 0),
        y=ag.obsm("pca", 1),
        color=ag.gene("gene_1", use_raw=False),
    )
    projected_values = adata.n_obs * 3

    with pytest.raises(ap.AnnplyrError, match=f"materialize {projected_values}"):
        ag.gganndata(
            adata,
            mapping,
            use_raw=False,
            max_matrix_values=projected_values - 1,
        )

    plot = ag.gganndata(
        adata,
        mapping,
        use_raw=False,
        max_matrix_values=projected_values,
    )
    assert list(plot.data.columns) == ["gene_1", "PC_1", "PC_2"]


def test_resolve_frame_handles_dataframe_backed_obsm_by_position():
    adata = _adata("dense")
    values = pd.DataFrame(
        np.arange(adata.n_obs * 3).reshape(adata.n_obs, 3),
        index=adata.obs_names,
        columns=["left", "middle", "right"],
    )
    adata.obsm["custom"] = values

    result = resolve_frame(
        adata,
        [ag.obsm("custom", 2), ag.obsm("custom", 0)],
        use_raw=False,
    )

    assert list(result.columns) == ["CUSTOM_3", "CUSTOM_1"]
    assert np.array_equal(result["CUSTOM_3"], values["right"])
    assert np.array_equal(result["CUSTOM_1"], values["left"])


def test_resolution_does_not_mutate_input_anndata():
    adata = _adata("csc")
    before_obs = adata.obs.copy(deep=True)
    before_var = adata.var.copy(deep=True)
    before_x = adata.X.copy()
    before_layer = adata.layers["double"].copy()
    before_raw = adata.raw.X.copy()
    before_obsm = adata.obsm["X_pca"].copy()
    before_uns = adata.uns["metadata"]["large"].copy()

    resolve_frame(
        adata,
        [
            ag.gene("gene_1", use_raw=False),
            ag.gene("gene_2", layer="double"),
            ag.gene("gene_3", use_raw=True),
            ag.obsm("pca", 4),
            "group",
        ],
        use_raw=False,
    )

    pd.testing.assert_frame_equal(adata.obs, before_obs)
    pd.testing.assert_frame_equal(adata.var, before_var)
    assert (adata.X != before_x).nnz == 0
    assert (adata.layers["double"] != before_layer).nnz == 0
    assert (adata.raw.X != before_raw).nnz == 0
    assert np.array_equal(adata.obsm["X_pca"], before_obsm)
    assert np.array_equal(adata.uns["metadata"]["large"], before_uns)


def test_source_validation_is_shared_by_accessors_and_helpers():
    adata = _adata("csr")
    with pytest.raises(KeyError, match="missing"):
        resolve_source(adata, "missing", None)
    with pytest.raises(ValueError, match="use_raw=True"):
        resolve_frame(
            adata,
            [ag.gene("gene_1", layer="double", use_raw=True)],
        )


def test_embedding_resolution_selects_requested_coordinates_in_order():
    adata = _adata("csr")
    result = resolve_frame(
        adata,
        [ag.obsm("pca", 7), ag.obsm("pca", 1), "group"],
        use_raw=False,
    )
    assert list(result.columns) == ["PC_8", "PC_2", "group"]
    assert np.allclose(result["PC_8"], adata.obsm["X_pca"][:, 7])
    assert np.allclose(result["PC_2"], adata.obsm["X_pca"][:, 1])


def test_resolve_frame_combines_x_layer_and_obs_sources():
    adata = _adata("csr")
    result = resolve_frame(
        adata,
        [
            ag.gene("gene_1", use_raw=False),
            ag.gene("gene_2", layer="double"),
            "group",
        ],
    )

    assert list(result.columns) == ["gene_1", "gene_2", "group"]
    assert np.allclose(result["gene_1"], adata.X[:, 1].toarray().ravel())
    assert np.allclose(result["gene_2"], adata.layers["double"][:, 2].toarray().ravel())


def test_aggregation_omits_empty_categories_without_mutating_order():
    adata = _adata("csc")
    adata.obs["group"] = adata.obs["group"].cat.add_categories(["empty"])
    categories_before = list(adata.obs["group"].cat.categories)

    result = aggregate_expression(adata, ["gene_1"], "group", use_raw=False)

    assert list(pd.unique(result["group"])) == ["a", "b", "c", "d"]
    assert list(adata.obs["group"].cat.categories) == categories_before


def test_aggregation_omits_rows_with_missing_group_keys():
    adata = _adata("csr")
    adata.obs.loc[adata.obs.index[:5], "group"] = np.nan

    result = aggregate_expression(adata, ["gene_1"], "group", use_raw=False)

    assert not result["group"].isna().any()
    retained = adata.obs["group"].notna().to_numpy()
    values = adata.X[retained, 1].toarray().ravel()
    groups = adata.obs.loc[retained, "group"]
    expected = pd.Series(values).groupby(groups.reset_index(drop=True), observed=True).mean()
    actual = result.set_index("group")["mean_expression"]
    assert np.allclose(actual.reindex(expected.index), expected)


def test_aggregation_temporary_columns_do_not_shadow_group_keys():
    adata = _adata("csr")
    adata.obs["__ggann_mean_0"] = adata.obs["group"]

    result = aggregate_expression(
        adata,
        ["gene_1"],
        "__ggann_mean_0",
        use_raw=False,
    )

    assert not result.empty
    assert "__ggann_mean_0" in result


def test_group_means_supports_obs_gene_name_collisions():
    adata = _adata("csc")
    adata.obs["gene_1"] = adata.obs["group"]

    result = group_means(adata, ["gene_1"], "gene_1", use_raw=False)

    assert list(result.columns) == ["gene_1"]
    assert len(result) == 4


def test_densify_frame_preserves_sparse_numeric_dtype():
    frame = pd.DataFrame.sparse.from_spmatrix(
        sparse.csr_matrix(np.eye(4, dtype=np.float32)), columns=list("abcd")
    )
    dense = densify_frame(frame)
    assert all(dtype == np.dtype("float32") for dtype in dense.dtypes)
    assert np.array_equal(dense.to_numpy(), np.eye(4, dtype=np.float32))


def test_highest_expression_sparse_and_dense_have_same_plot_data():
    dense = _adata("dense")
    csr = _adata("csr")
    dense_plot = ag.plot_highest_expr_genes(dense, n=6, use_raw=False)
    sparse_plot = ag.plot_highest_expr_genes(csr, n=6, use_raw=False)

    assert list(dense_plot.data["gene"].cat.categories) == list(
        sparse_plot.data["gene"].cat.categories
    )
    assert np.allclose(dense_plot.data["percent"], sparse_plot.data["percent"], equal_nan=True)
