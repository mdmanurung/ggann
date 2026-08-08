import pandas as pd
import pytest

from ggann import gene, obs, obsm
from ggann._resolve import embedding_coords, embedding_key, resolve_frame


def test_resolve_obs_and_gene(adata, group_key):
    df = resolve_frame(adata, [group_key, "CD3D"])
    assert list(df.columns) == [group_key, "CD3D"]
    assert len(df) == adata.n_obs
    assert pd.api.types.is_numeric_dtype(df["CD3D"])
    # gene column must be dense, not sparse
    assert not isinstance(df["CD3D"].dtype, pd.SparseDtype)


def test_resolve_embedding_names(adata):
    df = resolve_frame(adata, ["UMAP_1", "UMAP_2"])
    assert list(df.columns) == ["UMAP_1", "UMAP_2"]
    assert len(df) == adata.n_obs


def test_embedding_key_aliases(adata):
    assert embedding_key(adata, "umap") == "X_umap"
    assert embedding_key(adata, "X_umap") == "X_umap"
    assert embedding_key(adata, "UMAP") == "X_umap"


def test_embedding_coords_naming(adata):
    coords = embedding_coords(adata, "umap")
    assert list(coords.columns) == ["UMAP_1", "UMAP_2"]


@pytest.mark.parametrize("index", [True, 1.5, "1"])
def test_obsm_rejects_non_integer_indices(index):
    with pytest.raises(TypeError, match="integer"):
        obsm("umap", index)


def test_force_source_escapes(adata, group_key):
    # obs() and gene() force a source
    df_obs = resolve_frame(adata, [obs(group_key)])
    assert group_key in df_obs.columns
    df_gene = resolve_frame(adata, [gene("CD3D")])
    assert "CD3D" in df_gene.columns


def test_unknown_name_is_skipped(adata):
    # a name that matches nothing is silently dropped (assumed computed aesthetic)
    df = resolve_frame(adata, ["CD3D", "definitely_not_a_column"])
    assert "CD3D" in df.columns
    assert "definitely_not_a_column" not in df.columns


def test_missing_forced_gene_raises(adata):
    with pytest.raises(KeyError):
        resolve_frame(adata, [gene("NOT_A_GENE")])
