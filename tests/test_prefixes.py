"""Optional strict prefixes: obs: / gene: / obsm:, plus gene:NAME@layer."""

import numpy as np
import pytest
from plotnine import geom_point, ggplot

from ggann import aes, gganndata
from ggann._resolve import ObsmRef, Ref, parse_token, plain_name, resolve_frame


@pytest.fixture(scope="module")
def adata_layers(adata):
    a = adata.copy()
    a.layers["logcounts"] = a.X.copy() * 10.0
    return a


def test_parse_token_forms():
    assert parse_token("obs:phase") == Ref("phase", "obs")
    assert parse_token("gene:CD3D") == Ref("CD3D", "gene")
    assert parse_token("gene:CD3D@logcounts") == Ref("CD3D", "gene", layer="logcounts")
    assert parse_token("gene:CD3D@raw") == Ref("CD3D", "gene", use_raw=True)
    assert parse_token("gene:CD3D@X") == Ref("CD3D", "gene", use_raw=False)
    assert parse_token("obsm:umap[1]") == ObsmRef("umap", 1)
    # bare strings pass through untouched
    assert parse_token("CD3D") == "CD3D"
    assert parse_token("factor(louvain)") == "factor(louvain)"


def test_obsm_bad_token_raises():
    with pytest.raises(ValueError, match="obsm:umap"):
        parse_token("obsm:umap")  # missing [index]


def test_prefix_strict_obs_vs_gene(adata):
    # a colliding name resolves per the prefix, no ambiguity
    df_obs = resolve_frame(adata, ["obs:bulk_labels"])
    assert "bulk_labels" in df_obs.columns
    df_gene = resolve_frame(adata, ["gene:CD3D"])
    assert "CD3D" in df_gene.columns


def test_prefix_layer_selection(adata_layers):
    a = adata_layers
    df = resolve_frame(a, ["gene:CD3D@logcounts", "gene:CD8A@X"])
    xi, yi = a.var_names.get_loc("CD3D"), a.var_names.get_loc("CD8A")
    assert np.allclose(df["CD3D"], np.asarray(a.X[:, xi]).ravel() * 10, atol=1e-4)
    assert np.allclose(df["CD8A"], np.asarray(a.X[:, yi]).ravel(), atol=1e-4)


def test_obsm_prefix_coordinates(adata):
    df = resolve_frame(adata, ["obsm:umap[0]", "obsm:umap[1]"])
    assert list(df.columns) == ["UMAP_1", "UMAP_2"]
    # match the raw obsm array
    assert np.allclose(df["UMAP_1"], adata.obsm["X_umap"][:, 0], atol=1e-5)


def test_plain_name_rewrites(adata):
    assert plain_name(adata, "gene:CD3D@logcounts") == "CD3D"
    assert plain_name(adata, "obs:phase") == "phase"
    assert plain_name(adata, "obsm:umap[0]") == "UMAP_1"
    assert plain_name(adata, "CD3D") == "CD3D"  # bare unchanged


def test_gganndata_all_prefix_strings(adata_layers, tmp_path):
    p = (
        gganndata(
            adata_layers,
            aes("obsm:umap[0]", "obsm:umap[1]", color="gene:CD3D@logcounts", shape="obs:phase"),
        )
        + geom_point()
    )
    assert isinstance(p, ggplot)
    for col in ("UMAP_1", "UMAP_2", "CD3D", "phase"):
        assert col in p.data.columns
    out = tmp_path / "prefix.png"
    p.save(out, verbose=False, width=5, height=4, dpi=60)
    assert out.exists()


def test_auto_resolution_still_default(adata):
    # bare strings keep working exactly as before
    p = gganndata(adata, aes("UMAP_1", "UMAP_2", color="CD3D")) + geom_point()
    assert "CD3D" in p.data.columns
