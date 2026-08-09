"""Per-aesthetic layer/source selection via the gene() accessor."""

import numpy as np
import pytest
from plotnine import geom_point, ggplot

from ggann import aes, gene, gganndata, obs
from ggann._resolve import resolve_frame


@pytest.fixture(scope="module")
def adata_layers(adata):
    a = adata.copy()
    a.layers["counts"] = a.X.copy()
    a.layers["logcounts"] = a.X.copy() * 10.0
    return a


def test_gene_layer_picks_right_matrix(adata_layers):
    a = adata_layers
    df = resolve_frame(a, [gene("CD3D", layer="logcounts"), gene("CD8A", layer="counts")])
    xi, yi = a.var_names.get_loc("CD3D"), a.var_names.get_loc("CD8A")
    assert np.allclose(df["CD3D"], np.asarray(a.X[:, xi]).ravel() * 10, atol=1e-4)
    assert np.allclose(df["CD8A"], np.asarray(a.X[:, yi]).ravel(), atol=1e-4)


def test_bare_gene_inherits_plotwide_layer(adata_layers):
    a = adata_layers
    df = resolve_frame(a, ["CD3D"], layer="logcounts")
    xi = a.var_names.get_loc("CD3D")
    assert np.allclose(df["CD3D"], np.asarray(a.X[:, xi]).ravel() * 10, atol=1e-4)


def test_gene_use_raw_override(adata_layers):
    a = adata_layers
    # explicit use_raw on the accessor overrides a plot-wide layer for this gene
    df = resolve_frame(a, [gene("CD3D", use_raw=True), "CD8A"], layer="counts")
    xi = a.raw.var_names.get_loc("CD3D")
    assert np.allclose(df["CD3D"], np.asarray(a.raw.X[:, xi].todense()).ravel(), atol=1e-4)


def test_ref_in_aes_renders(adata_layers, tmp_path):
    a = adata_layers
    p = gganndata(a, aes("UMAP_1", "UMAP_2", color=gene("CD3D", layer="logcounts"))) + geom_point()
    assert isinstance(p, ggplot)
    assert "CD3D" in p.data.columns  # Ref resolved to a plain column
    out = tmp_path / "ref.png"
    p.save(out, verbose=False, width=5, height=4, dpi=60)
    assert out.exists()


def test_obs_accessor_in_aes_renders(adata_layers, tmp_path):
    p = gganndata(adata_layers, aes("UMAP_1", "UMAP_2", color=obs("phase"))) + geom_point()
    assert "phase" in p.data.columns
    out = tmp_path / "obsref.png"
    p.save(out, verbose=False, width=5, height=4, dpi=60)
    assert out.exists()
