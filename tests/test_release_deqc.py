"""Regressions for the 0.1.0 readiness rows owned by the ``deqc`` area.

Row IDs refer to plans/evaluation/remaining_todos.md. Each test name starts with
the row it closes.
"""

from __future__ import annotations

import inspect
import warnings

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotnine as p9
import pytest
import scanpy as sc
import scipy.stats as ss

import ggann as ag

GROUP = "bulk_labels"
KEY = "rank_wilcoxon"


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    plt.close("all")  # PyComplexHeatmap and marsilea draw on the current figure


@pytest.fixture(scope="module")
def de_adata(adata):
    """pbmc68k with a full-gene wilcoxon result under its own key."""
    a = adata.copy()
    sc.tl.rank_genes_groups(a, GROUP, method="wilcoxon", n_genes=a.raw.n_vars, key_added=KEY)
    return a


@pytest.fixture(scope="module")
def small_de(adata):
    a = adata.copy()
    sc.tl.rank_genes_groups(a, GROUP, method="wilcoxon", n_genes=20, key_added=KEY)
    return a


def _build(plot):
    assert isinstance(plot, p9.ggplot)
    plot._build()
    return plot


# --------------------------------------------------------------------------- M1


def test_C2_clustermap_numeric_annotations_of_any_width(adata):
    a = adata.copy()
    a.obs["f32"] = a.obs["n_counts"].astype("float32")
    a.obs["i32"] = np.arange(a.n_obs, dtype="int32")
    cm = ag.plot_clustermap(a, ["CD3D", "NKG7", "CST3"], None, annotations=["f32", "i32"])
    assert cm is not None


def test_C3_clustermap_checks_dimensions_against_enabled_clustering(adata, capsys):
    with pytest.raises(ValueError, match="row clustering needs at least 2 rows"):
        ag.plot_clustermap(adata, ["CD3D"], GROUP)
    # a single gene is fine once row clustering is off, even right after the failure
    assert ag.plot_clustermap(adata, ["CD3D"], GROUP, row_cluster=False) is not None
    # D10: PyComplexHeatmap progress text is off by default
    assert "Starting" not in capsys.readouterr().out


def test_C3_clustermap_single_group_needs_col_cluster_off(adata):
    one = adata[adata.obs[GROUP] == "CD56+ NK"].copy()
    one.obs[GROUP] = pd.Categorical(one.obs[GROUP].astype(str))
    with pytest.raises(ValueError, match="col clustering needs at least 2 cols"):
        ag.plot_clustermap(one, ["CD3D", "NKG7"], GROUP)
    assert ag.plot_clustermap(one, ["CD3D", "NKG7"], GROUP, col_cluster=False) is not None


@pytest.mark.parametrize("orientation", ["top", "left"])
def test_C4_dendrogram_draws_a_path_in_both_orientations(adata, orientation):
    plot = _build(ag.plot_dendrogram(adata, GROUP, orientation=orientation))
    assert type(plot.layers[0].geom).__name__ == "geom_path"


def _counts_adata():
    rng = np.random.default_rng(0)
    n, genes = 60, [f"g{i}" for i in range(8)]
    counts = rng.poisson(3, size=(n, len(genes))).astype(float)
    obs = pd.DataFrame(
        {
            "donor": np.repeat(["d1", "d2", "d3"], n // 3),
            "cell_type": np.tile(["A", "B"], n // 2),
            "condition": np.repeat(["x", "y", "x"], n // 3),
        },
        index=[f"c{i}" for i in range(n)],
    )
    a = ad.AnnData(counts[:, :5].copy(), obs=obs, var=pd.DataFrame(index=genes[:5]))
    a.raw = ad.AnnData(counts, obs=obs, var=pd.DataFrame(index=genes))
    return a, counts


def test_C6_pseudobulk_use_raw_with_different_raw_and_current_genes():
    a, counts = _counts_adata()
    pb = ag.pseudobulk(a, "donor", "cell_type", use_raw=True, min_cells=0)
    assert list(pb.var_names) == [f"g{i}" for i in range(8)]
    assert "condition" in pb.obs
    mask = (a.obs["donor"] == "d2") & (a.obs["cell_type"] == "B")
    row = (pb.obs["donor"] == "d2") & (pb.obs["cell_type"] == "B")
    np.testing.assert_allclose(np.asarray(pb[row].X).ravel(), counts[mask.to_numpy()].sum(axis=0))
    assert a.n_vars == 5  # input untouched


def test_C6_pseudobulk_use_raw_without_raw_is_a_clear_error():
    a, _ = _counts_adata()
    a.raw = None
    with pytest.raises(ValueError, match="adata.raw is None"):
        ag.pseudobulk(a, "donor", use_raw=True)


def _de_table(n=40):
    rng = np.random.default_rng(0)
    return pd.DataFrame(
        {
            "baseMean": rng.gamma(2.0, 50.0, n),
            "log2FoldChange": rng.normal(0, 2, n),
            "padj": rng.uniform(0, 1, n),
        },
        index=[f"G{i}" for i in range(n)],
    )


def test_C8_ma_warns_with_the_count_of_non_positive_means():
    de = _de_table()
    de.iloc[:3, 0] = [0.0, -1.0, np.nan]
    with pytest.warns(UserWarning, match=r"dropped 3 of 40 genes"):
        plot = ag.plot_ma(de)
    assert len(plot.data) == 37
    _build(plot)


def test_C8_ma_all_non_positive_raises():
    de = _de_table()
    de["baseMean"] = 0.0
    with pytest.warns(UserWarning, match="dropped 40 of 40"):
        with pytest.raises(ValueError, match="no gene has a positive"):
            ag.plot_ma(de)


@pytest.mark.parametrize("bad", [0, -3, 2.5, None, True])
@pytest.mark.parametrize("fn", ["plot_rank_genes_dotplot", "plot_rank_genes_matrixplot"])
def test_B6a_rank_genes_rejects_non_positive_n_genes(small_de, fn, bad):
    with pytest.raises(ValueError, match="n_genes"):
        getattr(ag, fn)(small_de, n_genes=bad, key=KEY)


def test_B6a_negative_n_genes_message_names_scanpy_limit(small_de):
    with pytest.raises(ValueError, match="bottom-ranked"):
        ag.plot_rank_genes_dotplot(small_de, n_genes=-3, key=KEY)


@pytest.mark.parametrize("bad", [-1, 2.5, True, "3"])
def test_B6b_ma_label_top_must_be_non_negative_integer(bad):
    with pytest.raises(ValueError, match="label_top"):
        ag.plot_ma(_de_table(), label_top=bad)


@pytest.mark.parametrize("bad", [0, -3, 2.5, None, True])
def test_B6c_variance_ratio_validates_n_pcs(adata, bad):
    with pytest.raises(ValueError, match="n_pcs"):
        ag.plot_variance_ratio(adata, n_pcs=bad)


@pytest.mark.parametrize("bad", [None, ["CD56+ NK"], ("CD56+ NK",), 3])
def test_B7_volcano_requires_one_str_group(small_de, bad):
    with pytest.raises(ValueError, match="one group name"):
        ag.plot_volcano(small_de, bad, key=KEY)


def test_B7_volcano_still_works_for_one_group(small_de):
    plot = ag.plot_volcano(small_de, "CD56+ NK", key=KEY)
    assert len(plot.data) == 20
    _build(plot)


# --------------------------------------------------------------------------- M2


def test_D2_scanpy_defaults():
    def default(fn, name):
        return inspect.signature(fn).parameters[name].default

    assert default(ag.plot_rank_genes_dotplot, "n_genes") == 10
    assert default(ag.plot_rank_genes_matrixplot, "n_genes") == 10
    assert default(ag.plot_rank_genes_matrixplot, "standard_scale") is None
    assert default(ag.plot_variance_ratio, "n_pcs") == 30
    assert default(ag.plot_variance_ratio, "log") is False


def test_D2_variance_ratio_defaults_draw_30_pcs_on_a_linear_axis(adata):
    a = adata.copy()
    a.uns["pca"]["variance_ratio"] = np.linspace(0.2, 0.001, 50)
    plot = ag.plot_variance_ratio(a)
    assert len(plot.data) == 30
    assert not any("log" in type(s).__name__ for s in plot.scales)


def test_D2_highest_expr_default_is_30_and_n_top_is_an_alias(adata):
    a = adata.copy()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        default = ag.plot_highest_expr_genes(a, use_raw=True)
        alias = ag.plot_highest_expr_genes(a, n_top=7, use_raw=True)
        same = ag.plot_highest_expr_genes(a, n=7, n_top=7, use_raw=True)
    assert default.data["gene"].nunique() == 30
    assert alias.data["gene"].nunique() == 7 == same.data["gene"].nunique()
    with pytest.raises(ValueError, match="conflict"):
        ag.plot_highest_expr_genes(a, n=5, n_top=7, use_raw=True)


def test_D2_rank_matrixplot_default_is_raw_group_means(small_de):
    plot = ag.plot_rank_genes_matrixplot(small_de, n_genes=3, key=KEY, groups=["CD19+ B"])
    manual = ag.plot_matrixplot(
        small_de, list(dict.fromkeys(plot.data["feature"].astype(str))), GROUP
    )
    np.testing.assert_allclose(plot.data["mean_expression"], manual.data["mean_expression"])


def test_D5_correlation_in_pca_space_matches_scanpy(adata):
    plot = ag.plot_correlation(adata, GROUP, use_rep="X_pca", n_pcs=10, cluster=False)
    got = plot.data.pivot(index="row", columns="col", values="corr").astype(float)
    pcs = pd.DataFrame(adata.obsm["X_pca"][:, :10], index=adata.obs_names)
    expected = pcs.groupby(adata.obs[GROUP].to_numpy(), observed=True).mean().T.corr()
    expected.columns = expected.columns.astype(str)
    expected.index = expected.index.astype(str)
    np.testing.assert_allclose(
        got.loc[expected.index, expected.columns].to_numpy(), expected.to_numpy(), atol=1e-6
    )


def test_D5_correlation_use_rep_conflicts(adata):
    with pytest.raises(ValueError, match="cannot be combined"):
        ag.plot_correlation(adata, GROUP, use_rep="X_pca", genes=["CD3D", "NKG7"])
    with pytest.raises(KeyError, match="use_rep"):
        ag.plot_correlation(adata, GROUP, use_rep="X_nope")
    with pytest.raises(ValueError, match="needs use_rep"):
        ag.plot_correlation(adata, GROUP, n_pcs=5, genes=["CD3D", "NKG7"])


def test_D9_ma_threshold_semantics_match_volcano():
    de = pd.DataFrame(
        {
            "baseMean": [10.0, 20.0, 30.0, 40.0, 50.0],
            "log2FoldChange": [2.0, -2.0, 0.5, 1.0, 3.0],
            "padj": [0.01, 0.01, 0.01, 0.2, 0.05],
        },
        index=list("abcde"),
    )
    plot = ag.plot_ma(de)
    sig = dict(zip(plot.data.index, plot.data["_significance"].astype(str)))
    # e: padj == 0.05 is not < 0.05; d: padj too large; c: |lfc| below 1
    assert sig == {
        "a": "up",
        "b": "down",
        "c": "not significant",
        "d": "not significant",
        "e": "not significant",
    }
    assert plot.data["significant"].sum() == 2
    assert ag.plot_ma(de, lfc=0).data["significant"].sum() == 3  # padj alone
    with pytest.raises(TypeError):
        ag.plot_ma(de, lfc="log2FoldChange")  # a column name is now lfc_col
    renamed = de.rename(columns={"log2FoldChange": "lfc"})
    assert len(ag.plot_ma(renamed, lfc_col="lfc").data) == 5
    _build(plot)


def test_D10_clustermap_documents_dendrogram_switch_and_passes_it_on(adata):
    assert "row_dendrogram=True" in inspect.getdoc(ag.plot_clustermap)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        cm = ag.plot_clustermap(
            adata, ["CD3D", "NKG7", "CST3"], GROUP, row_dendrogram=True, col_dendrogram=True
        )
    assert cm is not None
    assert not [w for w in caught if "set_bad" in str(w.message)]


# --------------------------------------------------------------------------- M3


def test_A5_groups_selects_marker_groups(small_de):
    names = small_de.uns[KEY]["names"]
    plot = ag.plot_rank_genes_matrixplot(small_de, n_genes=3, key=KEY, groups=["CD56+ NK"])
    assert list(plot.data["feature"].cat.categories) == list(names["CD56+ NK"][:3])
    with pytest.raises(KeyError, match="nope"):
        ag.plot_rank_genes_dotplot(small_de, key=KEY, groups=["nope"])


def test_A5_min_logfoldchange_filters_before_top_n(small_de):
    group = "CD56+ NK"
    table = ag.rank_genes_df(small_de, group=group, key=KEY)
    cutoff = float(table["logfoldchanges"].median())
    expected = list(table[table["logfoldchanges"] > cutoff]["names"][:4])
    plot = ag.plot_rank_genes_matrixplot(
        small_de, n_genes=4, key=KEY, groups=[group], min_logfoldchange=cutoff
    )
    assert list(plot.data["feature"].cat.categories) == expected


def test_A5_values_to_plot_uses_the_de_statistic_not_expression(de_adata):
    group = "CD56+ NK"
    plot = ag.plot_rank_genes_matrixplot(
        de_adata, n_genes=3, key=KEY, groups=[group], values_to_plot="scores"
    )
    scores = ag.rank_genes_df(de_adata, group=group, key=KEY).set_index("names")["scores"]
    shown = plot.data[plot.data[GROUP] == group]
    shown = shown.set_index(shown["feature"].astype(str))
    for gene, value in shown["mean_expression"].items():
        assert value == pytest.approx(scores[gene])
    # and it is not the expression mean
    means = ag.plot_matrixplot(de_adata, list(shown.index), GROUP).data
    assert not np.allclose(plot.data["mean_expression"], means["mean_expression"])
    log10 = ag.plot_rank_genes_dotplot(
        de_adata, n_genes=2, key=KEY, groups=[group], values_to_plot="log10_pvals_adj"
    )
    assert (log10.data["mean_expression"] >= 0).all()
    _build(log10)


def test_A5_values_to_plot_needs_all_genes_stored_and_no_standard_scale(small_de):
    with pytest.raises(ValueError, match="lacks statistics"):
        ag.plot_rank_genes_dotplot(small_de, n_genes=3, key=KEY, values_to_plot="scores")
    with pytest.raises(ValueError, match="values_to_plot"):
        ag.plot_rank_genes_dotplot(small_de, n_genes=3, key=KEY, values_to_plot="bogus")
    with pytest.raises(ValueError, match="standard_scale"):
        ag.plot_rank_genes_dotplot(
            small_de, n_genes=3, key=KEY, values_to_plot="scores", standard_scale="var"
        )


@pytest.mark.parametrize("fn", ["plot_rank_genes_dotplot", "plot_rank_genes_matrixplot"])
def test_A5_repeated_markers_are_one_labelled_column_per_occurrence(small_de, fn):
    a = small_de.copy()
    names = a.uns[KEY]["names"]
    shared = names["CD19+ B"][0]
    names["CD56+ NK"][0] = shared  # the same top marker in two groups
    plot = getattr(ag, fn)(a, n_genes=2, key=KEY, groups=["CD19+ B", "CD56+ NK"])
    cats = list(plot.data["feature"].cat.categories)
    assert cats == [
        f"{shared} (CD19+ B)",
        names["CD19+ B"][1],
        f"{shared} (CD56+ NK)",
        names["CD56+ NK"][1],
    ]
    one = plot.data[plot.data["feature"].astype(str).str.startswith(shared)]
    assert set(one["marker_group"]) == {"CD19+ B", "CD56+ NK"}
    # exact copies: one computation per gene
    a1, a2 = (one[one["marker_group"] == m].sort_values(GROUP) for m in ("CD19+ B", "CD56+ NK"))
    np.testing.assert_allclose(a1["mean_expression"].to_numpy(), a2["mean_expression"].to_numpy())
    _build(plot)


def test_A6_qc_violin_log_is_log1p_and_keeps_zeros(adata):
    a = adata.copy()
    a.obs.loc[a.obs_names[:5], "n_counts"] = 0.0
    raw = ag.plot_qc_violin(a, ["n_counts"])
    logged = ag.plot_qc_violin(a, ["n_counts"], log=True)
    np.testing.assert_allclose(logged.data["value"], np.log1p(raw.data["value"]))
    assert (logged.data["value"] == 0).sum() == 5
    assert logged.labels.y == "log1p(value)"
    a.obs["neg"] = -1.0
    with pytest.raises(ValueError, match="non-negative"):
        ag.plot_qc_violin(a, ["neg"], log=True)


def test_A6_qc_scatter_alpha_and_log(adata):
    plot = ag.plot_qc_scatter(adata, "n_counts", "n_genes", alpha=0.2, log=True)
    assert plot.layers[-1].geom.aes_params["alpha"] == 0.2
    assert [type(s).__name__ for s in plot.scales] == ["scale_x_log10", "scale_y_log10"]
    with pytest.raises(ValueError, match="alpha"):
        ag.plot_qc_scatter(adata, "n_counts", "n_genes", alpha=0)
    a = adata.copy()
    a.obs.loc[a.obs_names[:4], "n_genes"] = 0
    with pytest.raises(ValueError, match="drop 4 cells"):
        ag.plot_qc_scatter(a, "n_counts", "n_genes", log=True)


def test_A6_highest_expr_gene_symbols_and_log(adata):
    a = adata.copy()
    a.raw.var["symbol"] = [f"S_{g}" for g in a.raw.var_names]
    a.raw.var["dup"] = "same"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        plain = ag.plot_highest_expr_genes(a, n=5, use_raw=True)
        named = ag.plot_highest_expr_genes(a, n=5, use_raw=True, gene_symbols="symbol")
        logged = ag.plot_highest_expr_genes(a, n=5, use_raw=True, log=True)
    assert [f"S_{g}" for g in plain.data["gene"].cat.categories] == list(
        named.data["gene"].cat.categories
    )
    np.testing.assert_allclose(logged.data["percent"], np.log1p(plain.data["percent"]))
    with pytest.raises(ValueError, match="duplicate"):
        ag.plot_highest_expr_genes(a, n=5, use_raw=True, gene_symbols="dup")
    with pytest.raises(KeyError, match="gene_symbols"):
        ag.plot_highest_expr_genes(a, n=5, use_raw=True, gene_symbols="missing")


# ------------------------------------------------------------------------ minor


def test_Q5_upset_sorts_largest_first_and_rejects_empty_sets():
    sets = {"A": list("abcdefgh"), "B": list("defghijk"), "C": list("hxyz")}

    def bar_heights(**kw):
        up = ag.plot_upset(sets, **kw)
        for axis in up.figure.axes:
            bars = sorted(
                (p.get_x(), int(p.get_height()))
                for p in axis.patches
                if p.get_height() > 0 and p.get_width() < 1.5
            )
            if bars:
                return [h for _, h in bars]
        raise AssertionError("no bars")

    default = bar_heights()
    assert default == sorted(default, reverse=True)
    assert bar_heights(ascending=True) == sorted(default)
    with pytest.raises(ValueError, match="every set is empty"):
        ag.plot_upset({"A": [], "B": []})
    with pytest.raises(ValueError, match="sort_subsets"):
        ag.plot_upset(sets, sort_subsets="size")


def test_Q6_pseudobulk_filter_is_reported_and_missing_labels_are_clear():
    a, _ = _counts_adata()
    with pytest.warns(UserWarning, match=r"dropped \d+ of \d+ pseudobulk profiles"):
        ag.pseudobulk(a, "donor", "cell_type", use_raw=True, min_cells=15)
    with pytest.raises(KeyError, match="'zzz' is not a column"):
        ag.pseudobulk(a, "zzz")
    with pytest.raises(KeyError, match="'zzz' is not a column"):
        ag.pseudobulk(a, "donor", "zzz")
    with pytest.raises(ValueError, match="min_cells"):
        ag.pseudobulk(a, "donor", min_cells=-1)
    nan = a.copy()
    nan.obs["donor"] = nan.obs["donor"].astype(object)
    nan.obs.iloc[:3, nan.obs.columns.get_loc("donor")] = np.nan
    with pytest.warns(UserWarning, match="dropped 3 of 60 cells with missing 'donor'"):
        ag.pseudobulk(nan, "donor", "cell_type", min_cells=0)


def test_Q9_dendrogram_single_group_unused_category_and_dtype(adata):
    one = adata.copy()
    one.obs["one"] = pd.Categorical(["x"] * one.n_obs)
    with pytest.raises(ValueError, match="needs at least 2"):
        ag.plot_dendrogram(one, "one")
    with pytest.raises(ValueError, match="categorical"):
        ag.plot_dendrogram(adata, "n_counts")
    with pytest.raises(KeyError, match="not a column"):
        ag.plot_dendrogram(adata, "nope")
    ghost = adata.copy()
    cats = [*ghost.obs[GROUP].cat.categories, "ghost"]
    ghost.obs[GROUP] = pd.Categorical(ghost.obs[GROUP], categories=cats)
    plot = _build(ag.plot_dendrogram(ghost, GROUP))
    assert len(plot.data) > 0
    assert list(ghost.obs[GROUP].cat.categories)[-1] == "ghost"  # input untouched


def test_Q10_ma_rejects_non_dataframes_and_missing_label_column():
    with pytest.raises(TypeError, match="DataFrame"):
        ag.plot_ma({"baseMean": [1.0]})
    with pytest.raises(KeyError, match="nope"):
        ag.plot_ma(_de_table(), label="nope")


def test_Q10_highest_expr_scaled_x_hint_names_raw(adata):
    with pytest.warns(UserWarning, match="adata.raw is available: pass use_raw=True"):
        ag.plot_highest_expr_genes(adata, n=3)


def test_Q17_qc_violin_rejects_non_numeric_metrics(adata):
    with pytest.raises(ValueError, match="numeric"):
        ag.plot_qc_violin(adata, [GROUP])


def test_Q17_volcano_rejects_colour_keyword(small_de):
    with pytest.raises(ValueError, match="direction"):
        ag.plot_volcano(small_de, "CD56+ NK", key=KEY, color="red")


def test_Q14_correlation_does_not_warn_for_a_small_panel(adata):
    a = adata[:, :50].copy()
    a.raw = None
    a.var.pop("highly_variable") if "highly_variable" in a.var else None
    with warnings.catch_warnings():
        warnings.simplefilter("error", UserWarning)
        ag.plot_correlation(a, GROUP, cluster=False)


# --------------------------------------------------------------------- VERIFY-4


@pytest.fixture(scope="module")
def stat_frame():
    rng = np.random.default_rng(1)
    x = rng.normal(size=60)
    return pd.DataFrame(
        {"x": x, "y": 2 * x + rng.normal(size=60), "g": np.repeat(["a", "b", "c"], 20)}
    )


def _layer_data(frame, mapping, layer):
    plot = p9.ggplot(frame, mapping) + p9.geom_point() + layer
    plot._build()
    return plot.layers[-1].data.iloc[0]


@pytest.mark.parametrize(
    ("method", "scipy_fn"), [("pearson", ss.pearsonr), ("spearman", ss.spearmanr)]
)
def test_VERIFY4_stat_cor_matches_scipy(stat_frame, method, scipy_fn):
    got = _layer_data(stat_frame, p9.aes("x", "y"), ag.stat_cor(method=method))
    expected = scipy_fn(stat_frame["x"], stat_frame["y"])
    assert got["r"] == pytest.approx(expected[0], abs=1e-9)
    assert got["p"] == pytest.approx(expected[1], rel=1e-6, abs=1e-12)


def test_VERIFY4_stat_regline_equation_matches_linregress(stat_frame):
    got = _layer_data(stat_frame, p9.aes("x", "y"), ag.stat_regline_equation())
    fit = ss.linregress(stat_frame["x"], stat_frame["y"])
    assert got["rr"] == pytest.approx(fit.rvalue**2, abs=1e-9)
    n = len(stat_frame)
    assert got["adj_rr"] == pytest.approx(1 - (1 - fit.rvalue**2) * (n - 1) / (n - 2), abs=1e-9)
    sign = "-" if fit.intercept < 0 else "+"
    assert got["eq"] == f"y = {fit.slope:.2g}x {sign} {abs(fit.intercept):.2g}"


def test_VERIFY4_stat_compare_means_matches_scipy(stat_frame):
    two = stat_frame[stat_frame["g"] != "c"]
    a, b = (two.loc[two["g"] == k, "y"] for k in ("a", "b"))
    mapping = p9.aes("g", "y")
    wilcoxon = _layer_data(two, mapping, ag.stat_compare_means(method="wilcox.test"))
    # R's wilcox.test is exact below 50 observations; scipy's asymptotic p with
    # continuity correction is what the layer reports, and "auto" selects it here.
    assert wilcoxon["p"] == pytest.approx(ss.mannwhitneyu(a, b).pvalue, rel=1e-6)
    welch = _layer_data(two, mapping, ag.stat_compare_means(method="t.test"))
    assert welch["p"] == pytest.approx(ss.ttest_ind(a, b, equal_var=False).pvalue, rel=1e-6)
    kruskal = _layer_data(stat_frame, mapping, ag.stat_compare_means(method="kruskal.test"))
    expected = ss.kruskal(*[stat_frame.loc[stat_frame["g"] == k, "y"] for k in "abc"])
    assert kruskal["p"] == pytest.approx(expected.pvalue, rel=1e-6)
