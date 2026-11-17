"""The aggregation must match scanpy's dotplot, not merely be internally consistent.

``sc.pl.dotplot`` aggregates ``adata.raw`` by default, so ggann's default
(``use_raw=None`` -> raw) must reproduce its ``dot_color_df`` (mean expression)
and ``dot_size_df`` (fraction expressing).
"""

import numpy as np

from ggann._aggregate import aggregate_expression
from ggann._expression import resolve_source


def test_source_defaults_to_raw(adata):
    kind, lyr = resolve_source(adata, layer=None, use_raw=None)
    assert kind == "raw"
    assert lyr is None


def test_matches_scanpy_dotplot(adata, markers, group_key):
    import scanpy as sc

    dp = sc.pl.dotplot(adata, markers, groupby=group_key, return_fig=True)
    sc_mean = dp.dot_color_df[markers]
    sc_frac = dp.dot_size_df[markers]

    agg = aggregate_expression(adata, markers, group_key)  # default -> raw
    mean = agg.pivot(index=group_key, columns="feature", values="mean_expression")
    frac = agg.pivot(index=group_key, columns="feature", values="fraction")

    mean = mean.reindex(sc_mean.index)[markers]
    frac = frac.reindex(sc_frac.index)[markers]

    assert np.allclose(mean.values, sc_mean.values, atol=1e-5)
    assert np.allclose(frac.values, sc_frac.values, atol=1e-6)


def test_standard_scale_var_bounded(adata, markers, group_key):
    agg = aggregate_expression(adata, markers, group_key, standard_scale="var")
    vals = agg["mean_expression"].to_numpy()
    assert vals.min() >= -1e-9
    assert vals.max() <= 1 + 1e-9


def test_feature_order_preserved(adata, markers, group_key):
    agg = aggregate_expression(adata, markers, group_key)
    assert list(agg["feature"].cat.categories) == list(markers)


def test_group_order_follows_declared_categories_not_row_order(adata, group_key):
    """The same cells in another row order must give the same group order."""
    import numpy as np

    import ggann as ag

    genes = list(adata.raw.var_names[:6])
    shuffled = adata[np.random.default_rng(0).permutation(adata.n_obs)].copy()
    declared = [c for c in adata.obs[group_key].cat.categories if c in set(adata.obs[group_key])]

    corr = ag.plot_correlation(shuffled, group_key, genes=genes, cluster=False)
    assert list(dict.fromkeys(corr.data["row"].astype(str))) == declared

    means = aggregate_expression(shuffled, genes, group_key, use_raw=True)
    assert list(dict.fromkeys(means[group_key].astype(str))) == declared
