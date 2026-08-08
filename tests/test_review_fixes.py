"""Regression tests for defects surfaced by the adversarial review."""

import numpy as np
import pandas as pd
import pytest

import ggann as ag
from ggann._resolve import resolve_frame, resolve_source


def _no_raw(adata):
    a = adata.copy()
    a.raw = None
    return a


def test_use_raw_true_without_raw_raises_everywhere(adata):
    a = _no_raw(adata)
    # grammar path (previously silently fell back to .X)
    with pytest.raises(ValueError, match="adata.raw is None"):
        resolve_frame(a, ["CD3D"], use_raw=True)
    # aggregation path
    with pytest.raises(ValueError, match="adata.raw is None"):
        resolve_source(a, None, True)


def test_layer_plus_use_raw_conflicts(adata):
    with pytest.raises(ValueError, match="use_raw=True and a layer"):
        resolve_source(adata, "somelayer", True)


def test_gene_only_in_X_warns_when_reading_raw(adata):
    # a gene present in .X but not in .raw should warn (not silently vanish)
    a = adata.copy()
    gene_x = a.var_names[0]  # in .X
    raw_ad = a.raw.to_adata()
    a.raw = raw_ad[:, [g for g in raw_ad.var_names if g != gene_x]].copy()
    assert gene_x in a.var_names and gene_x not in set(a.raw.var_names)
    with pytest.warns(UserWarning, match="active expression source"):
        df = resolve_frame(a, [gene_x])  # defaults to raw
    assert gene_x not in df.columns


def test_group_order_follows_obs_categorical(adata, markers, group_key):
    # integer-named (Leiden-style) categories must not sort lexicographically
    a = adata.copy()
    a.obs["clust"] = pd.Categorical(
        (a.obs[group_key].cat.codes % 12).astype(str),
        categories=[str(i) for i in range(12)],
        ordered=True,
    )
    p = ag.plot_dotplot(a, markers, "clust")
    order = list(p.data["clust"].cat.categories)
    present = [str(i) for i in range(12) if str(i) in order]
    assert order == present  # numeric order, not '0','1','10','11','2'...


def test_group_by_is_unordered_categorical(adata, markers, group_key):
    # nominal identities -> unordered categorical -> qualitative palette on violin fill
    p = ag.plot_violin(adata, markers, group_key)
    assert p.data[group_key].cat.ordered is False


def test_duplicate_categories_order_is_deduped(adata, markers, group_key):
    cats = list(adata.obs[group_key].cat.categories)
    p = ag.plot_dotplot(adata, markers, group_key, categories_order=cats + cats)
    assert len(set(p.data[group_key].cat.categories)) == len(
        p.data[group_key].cat.categories
    )


def test_array_categories_order_is_supported(adata, markers, group_key):
    categories = np.asarray(adata.obs[group_key].cat.categories)
    plots = [
        ag.plot_stacked_violin(
            adata, markers[:1], group_key, categories_order=categories
        ),
        ag.plot_tracksplot(adata, markers[:1], group_key, categories_order=categories),
        ag.plot_proportions(
            adata,
            group_key,
            split_by="phase",
            categories_order=categories,
        ),
    ]
    for plot in plots:
        assert list(plot.data[group_key].cat.categories) == list(categories)


@pytest.mark.parametrize(
    "builder",
    [
        lambda ad, genes, group: ag.plot_stacked_violin(
            ad, genes[:1], group, categories_order=[]
        ),
        lambda ad, genes, group: ag.plot_tracksplot(
            ad, genes[:1], group, categories_order=[]
        ),
        lambda ad, genes, group: ag.plot_proportions(
            ad, group, split_by="phase", categories_order=[]
        ),
    ],
)
def test_empty_categories_order_is_not_treated_as_none(
    adata, markers, group_key, builder
):
    with pytest.raises(ValueError, match="missing groups"):
        builder(adata, markers, group_key)


def test_missing_group_values_are_preserved_for_plotnine_to_drop(
    adata, markers, group_key
):
    ad = adata.copy()
    ad.obs[group_key] = ad.obs[group_key].astype(object)
    ad.obs.iloc[0, ad.obs.columns.get_loc(group_key)] = None

    plot = ag.plot_violin(ad, markers[:1], group_key)

    assert plot.data[group_key].isna().sum() == 1


def test_partial_categories_order_raises(adata, markers, group_key):
    cats = list(adata.obs[group_key].cat.categories)[:2]
    with pytest.raises(ValueError, match="missing groups"):
        ag.plot_dotplot(adata, markers, group_key, categories_order=cats)


def test_matrixplot_default_matches_scanpy_scale(adata, markers, group_key):
    import scanpy as sc

    # default standard_scale=None -> raw group means, comparable to sc.pl.matrixplot
    p = ag.plot_matrixplot(adata, markers, group_key)
    got = p.data.pivot_table(
        index=group_key, columns="feature", values="mean_expression", observed=True
    )
    dp = sc.pl.dotplot(adata, markers, groupby=group_key, return_fig=True)
    sc_mean = dp.dot_color_df[markers].reindex(got.index)[markers]
    assert np.allclose(got[markers].values, sc_mean.values, atol=1e-4)


def test_embedding_requires_2d(adata):
    a = adata.copy()
    a.obsm["X_1d"] = a.obsm["X_umap"][:, :1]
    with pytest.raises(ValueError, match="at least 2"):
        ag.plot_embedding(a, basis="1d")


def test_clustermap_standard_scale_and_zscore_mutually_exclusive(
    adata, markers, group_key
):
    pytest.importorskip("PyComplexHeatmap")
    with pytest.raises(ValueError, match="mutually exclusive"):
        ag.plot_clustermap(
            adata, markers, group_by=group_key, standard_scale="var", z_score=0
        )
