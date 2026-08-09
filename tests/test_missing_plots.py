"""The missing-plots backlog surfaced by the scanpy comparison."""

from __future__ import annotations

import numpy as np
import pandas as pd
import plotnine as p9
import pytest

import ggann as ag


def test_variance_ratio_builds(adata):
    plot = ag.plot_variance_ratio(adata, n_pcs=30)
    assert isinstance(plot, p9.ggplot)
    assert len(plot.data) == 30  # capped at n_pcs
    plot._build()


def test_variance_ratio_missing_pca_raises(adata):
    ad = adata.copy()
    ad.uns.pop("pca", None)
    with pytest.raises(KeyError):
        ag.plot_variance_ratio(ad)


def test_embedding_density_grouped_and_pooled(adata, group_key):
    grouped = ag.plot_embedding_density(adata, "umap", group_key)
    pooled = ag.plot_embedding_density(adata, "umap")
    for plot in (grouped, pooled):
        assert isinstance(plot, p9.ggplot)
        # density is min-max scaled into [0, 1]
        assert plot.data["density"].between(0, 1).all()
        plot._build()


def test_embedding_density_assigns_each_cell_its_own_density(group_key):
    """Each cell's density must come from its OWN group, in cell order. A
    degenerate 2-cell group B has density exactly 0 (KDE needs >2 points); if the
    grouped result is assigned positionally instead of by index, those zeros land
    on the wrong (interleaved) rows and the B cells pick up nonzero A densities."""
    import anndata as ad_mod
    import numpy as np
    import pandas as pd

    rng = np.random.RandomState(0)
    n = 42
    groups = np.array(["A"] * n, dtype=object)
    groups[1] = groups[3] = "B"  # two B cells, interleaved among the A cells
    coords = rng.normal(0, 0.05, (n, 2))  # tight A blob
    coords[[1, 3]] = [[10.0, 10.0], [10.1, 9.9]]  # B lives far away
    coords[20] = [5.0, 5.0]  # an A outlier, for the secondary check
    adata = ad_mod.AnnData(
        X=rng.normal(size=(n, 3)),
        obs=pd.DataFrame({group_key: pd.Categorical(groups)}, index=[f"c{i}" for i in range(n)]),
    )
    adata.obsm["X_umap"] = coords
    d = ag.plot_embedding_density(adata, "umap", group_key).data

    # the degenerate B group must be exactly 0 everywhere (the sharp discriminator)
    assert (d.loc[d[group_key] == "B", "density"] == 0).all()
    # and within A, the far outlier is less dense than the blob median
    a = d[d[group_key] == "A"]
    xcol = next(c for c in d.columns if str(c).startswith("UMAP"))
    outlier = a.loc[a[xcol].idxmax()]
    assert outlier["density"] < a.drop(index=outlier.name)["density"].median()


def test_heatmap_builds_and_scales(adata, markers, group_key):
    plain = ag.plot_heatmap(adata, markers, group_key, use_raw=True)
    scaled = ag.plot_heatmap(adata, markers, group_key, use_raw=True, standard_scale="var")
    assert isinstance(plain, p9.ggplot)
    plain._build()
    # each gene scaled independently into [0, 1]
    assert scaled.data["value"].between(0, 1).all()
    scaled._build()


def test_dendrogram_builds_and_autocomputes(adata, group_key):
    ad = adata.copy()
    ad.uns.pop(f"dendrogram_{group_key}", None)
    plot = ag.plot_dendrogram(ad, group_key)
    assert isinstance(plot, p9.ggplot)
    assert f"dendrogram_{group_key}" not in ad.uns
    plot._build()
    ag.plot_dendrogram(ad, group_key, orientation="left")._build()


def test_dendrogram_bad_orientation(adata, group_key):
    with pytest.raises(ValueError):
        ag.plot_dendrogram(adata, group_key, orientation="sideways")


def test_dendrogram_custom_key_autocomputes_into_that_key(adata, group_key):
    ad = adata.copy()
    ad.uns.pop(f"dendrogram_{group_key}", None)
    plot = ag.plot_dendrogram(ad, group_key, key="my_dendro")
    assert "my_dendro" not in ad.uns
    plot._build()


def test_sina_builds_multi_and_single(adata, markers, group_key):
    multi = ag.plot_sina(adata, markers, group_key, use_raw=True, downsample=100)
    single = ag.plot_sina(adata, markers[:1], group_key, use_raw=True, violin=False)
    for plot in (multi, single):
        assert isinstance(plot, p9.ggplot)
        plot._build()


def _fake_de(n=400, seed=0):
    rng = np.random.RandomState(seed)
    return pd.DataFrame(
        {
            "baseMean": rng.gamma(2.0, 50.0, n),
            "log2FoldChange": rng.normal(0, 2, n),
            "padj": rng.uniform(0, 1, n),
        },
        index=[f"G{i}" for i in range(n)],
    )


def test_ma_builds_and_flags_significant():
    de = _fake_de()
    plot = ag.plot_ma(de, label_top=5)
    assert isinstance(plot, p9.ggplot)
    expected = ((de["baseMean"] > 0) & (de["padj"] < 0.05)).sum()
    assert plot.data["significant"].sum() == expected
    plot._build()


def test_ma_missing_columns_raises():
    with pytest.raises(KeyError):
        ag.plot_ma(pd.DataFrame({"baseMean": [1.0], "log2FoldChange": [0.5]}))
