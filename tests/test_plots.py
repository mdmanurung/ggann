import pytest
from plotnine import geom_point, ggplot

import ggann as ag
from ggann import aes, gganndata


def _renders(plot, tmp_path, name):
    """A ggplot is only 'valid' if it can actually draw; save to a temp PNG."""
    out = tmp_path / f"{name}.png"
    plot.save(out, verbose=False, width=5, height=4, dpi=60)
    assert out.exists() and out.stat().st_size > 0


def test_gganndata_returns_ggplot(adata):
    p = gganndata(adata, aes("UMAP_1", "UMAP_2", color="bulk_labels")) + geom_point()
    assert isinstance(p, ggplot)


def test_gganndata_resolves_gene(adata, tmp_path):
    p = gganndata(adata, aes("UMAP_1", "UMAP_2", color="CD3D")) + geom_point()
    assert "CD3D" in p.data.columns
    _renders(p, tmp_path, "gg_gene")


def test_plot_embedding_density(adata, tmp_path):
    p = ag.plot_embedding(adata, basis="umap")  # no color -> pointdensity
    assert isinstance(p, ggplot)
    _renders(p, tmp_path, "emb_density")


def test_plot_embedding_categorical(adata, tmp_path):
    p = ag.plot_embedding(adata, basis="umap", color="bulk_labels")
    assert isinstance(p, ggplot)
    _renders(p, tmp_path, "emb_cat")


def test_plot_embedding_gene(adata, tmp_path):
    p = ag.plot_embedding(adata, basis="umap", color="CD3D")
    assert isinstance(p, ggplot)
    _renders(p, tmp_path, "emb_gene")


def test_plot_dotplot(adata, markers, group_key, tmp_path):
    p = ag.plot_dotplot(adata, markers, group_key)
    assert isinstance(p, ggplot)
    _renders(p, tmp_path, "dotplot")


def test_plot_matrixplot(adata, markers, group_key, tmp_path):
    p = ag.plot_matrixplot(adata, markers, group_key)
    assert isinstance(p, ggplot)
    _renders(p, tmp_path, "matrixplot")


def test_plot_violin(adata, markers, group_key, tmp_path):
    p = ag.plot_violin(adata, markers, group_key)
    assert isinstance(p, ggplot)
    _renders(p, tmp_path, "violin")


def test_plot_clustermap(adata, markers, group_key):
    pch = pytest.importorskip("PyComplexHeatmap")
    import matplotlib.pyplot as plt

    plt.figure()
    cm = ag.plot_clustermap(adata, markers, group_by=group_key)
    assert isinstance(cm, pch.ClusterMapPlotter)
    plt.close("all")
