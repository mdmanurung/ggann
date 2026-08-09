"""Contracts for the explicit, composable Matplotlib acceleration path."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData, read_h5ad
from matplotlib.collections import PathCollection, QuadMesh
from pandas.testing import assert_frame_equal
from PIL import Image
from plotnine import ggplot, theme
from scipy import sparse

import ggann as ag
from benchmarks.compare_scanpy import _adata_digest


def _plots(adata, markers, group_key):
    return {
        "embedding": lambda backend: ag.plot_embedding(
            adata,
            "umap",
            color=group_key,
            pointdensity=False,
            backend=backend,
        ),
        "dotplot": lambda backend: ag.plot_dotplot(
            adata,
            markers,
            group_key,
            backend=backend,
        ),
        "matrixplot": lambda backend: ag.plot_matrixplot(
            adata,
            markers,
            group_key,
            backend=backend,
        ),
    }


def _source_adata(matrix_format: str) -> AnnData:
    values = np.array(
        [
            [0.0, 1.0, np.nan, 4.0],
            [2.0, 0.0, 3.0, 0.0],
            [0.0, 5.0, 0.0, 1.0],
            [4.0, 0.0, 2.0, np.nan],
            [1.0, 2.0, 0.0, 3.0],
            [0.0, 0.0, 1.0, 0.0],
        ],
        dtype=np.float32,
    )
    matrix = {
        "dense": values,
        "csr_matrix": sparse.csr_matrix(values),
        "csc_matrix": sparse.csc_matrix(values),
        "csr_array": sparse.csr_array(values),
        "csc_array": sparse.csc_array(values),
    }[matrix_format]
    result = AnnData(
        matrix,
        obs=pd.DataFrame(
            {"group": pd.Categorical(["b", "a", None, "b", "a", "b"])},
            index=["dup", "dup", "c2", "c3", "c4", "c5"],
        ),
        var=pd.DataFrame(index=["g0", "g1", "g2", "g3"]),
    )
    result.layers["double"] = matrix * 2
    result.raw = result.copy()
    result.obsm["X_umap"] = np.arange(result.n_obs * 2, dtype=np.float32).reshape(-1, 2)
    return result


@pytest.mark.parametrize("name", ["embedding", "dotplot", "matrixplot"])
def test_fast_backend_preserves_grammar_payload(adata, markers, group_key, name):
    default = _plots(adata, markers, group_key)[name]("plotnine")
    direct = _plots(adata, markers, group_key)[name]("matplotlib")

    assert isinstance(direct, ggplot)
    assert isinstance(direct, ag.MatplotlibGGPlot)
    assert direct.fast_path_active
    assert list(direct.mapping) == list(default.mapping)
    assert direct.labels == default.labels
    assert_frame_equal(
        direct.data.reset_index(drop=True),
        default.data.reset_index(drop=True),
        check_dtype=True,
        check_categorical=True,
        rtol=1e-7,
        atol=1e-8,
    )


@pytest.mark.parametrize("name", ["embedding", "dotplot", "matrixplot"])
def test_fast_backend_draws_png_and_vector_svg(
    adata,
    markers,
    group_key,
    tmp_path: Path,
    name,
):
    plot = _plots(adata, markers, group_key)[name]("matplotlib")
    figure = plot.draw(show=False)
    png = tmp_path / f"{name}.png"
    svg = tmp_path / f"{name}.svg"
    figure.savefig(png, format="png", dpi=80)
    figure.savefig(svg, format="svg")

    image = np.asarray(Image.open(png).convert("RGB"), dtype=np.float32) / 255
    foreground = np.mean(np.any(image < 0.98, axis=2))
    assert image.shape[:2] == (384, 512)
    assert float(image.std()) > 0.05
    assert 0.005 < foreground < 0.8
    assert "<image" not in svg.read_text()

    main = figure.axes[0]
    if name in {"embedding", "dotplot"}:
        assert sum(isinstance(item, PathCollection) for item in main.collections) == 1
    else:
        assert sum(isinstance(item, QuadMesh) for item in main.collections) == 1
        assert all(not item.get_rasterized() for item in main.collections)


def test_grammar_composition_disables_only_the_new_fast_plot(adata, group_key):
    direct = ag.plot_embedding(
        adata,
        color=group_key,
        pointdensity=False,
        backend="matplotlib",
    )
    composed = direct + theme(figure_size=(5, 4))

    assert direct.fast_path_active
    assert isinstance(composed, ag.MatplotlibGGPlot)
    assert not composed.fast_path_active
    figure = composed.draw(show=False)
    assert tuple(np.round(figure.get_size_inches(), 6)) == (5.0, 4.0)


def test_fast_backend_continuous_color_and_input_immutability(adata):
    before = _adata_digest(adata)
    plot = ag.plot_embedding(
        adata,
        color="n_genes",
        pointdensity=False,
        backend="matplotlib",
    )
    figure = plot.draw(show=False)
    assert len(figure.axes) == 2
    assert pd.api.types.is_numeric_dtype(plot.data["n_genes"])
    assert _adata_digest(adata) == before


@pytest.mark.parametrize(
    "matrix_format",
    ["dense", "csr_matrix", "csc_matrix", "csr_array", "csc_array"],
)
@pytest.mark.parametrize(
    "source_kwargs",
    [{"use_raw": False}, {"use_raw": True}, {"layer": "double"}],
)
def test_native_grouped_payload_matches_grammar_for_sources_and_sparse_types(
    matrix_format,
    source_kwargs,
):
    with pytest.warns(UserWarning, match="Observation names are not unique"):
        adata = _source_adata(matrix_format)
    for constructor in (ag.plot_dotplot, ag.plot_matrixplot):
        default = constructor(adata, ["g3", "g1", "g2"], "group", **source_kwargs)
        direct = constructor(
            adata,
            ["g3", "g1", "g2"],
            "group",
            backend="matplotlib",
            **source_kwargs,
        )
        assert_frame_equal(
            direct.data.reset_index(drop=True),
            default.data.reset_index(drop=True),
            check_dtype=True,
            check_categorical=True,
            rtol=1e-7,
            atol=1e-8,
        )


@pytest.mark.parametrize("matrix_format", ["dense", "csr_matrix", "csc_matrix"])
def test_native_backend_matches_backed_and_reordered_view(tmp_path, matrix_format):
    with pytest.warns(UserWarning, match="Observation names are not unique"):
        adata = _source_adata(matrix_format)
    expected = ag.plot_dotplot(
        adata[[5, 1, 3, 0]],
        ["g3", "g0"],
        "group",
        use_raw=False,
        backend="matplotlib",
    ).data
    path = tmp_path / f"direct-{matrix_format}.h5ad"
    adata.write_h5ad(path)
    backed = read_h5ad(path, backed="r")
    try:
        actual = ag.plot_dotplot(
            backed[[5, 1, 3, 0]],
            ["g3", "g0"],
            "group",
            use_raw=False,
            backend="matplotlib",
        ).data
        assert backed.isbacked
    finally:
        backed.file.close()
    assert_frame_equal(actual, expected, check_dtype=False, rtol=1e-7, atol=1e-8)


@pytest.mark.parametrize(
    ("call", "match"),
    [
        (lambda ad, genes, group: ag.plot_embedding(ad, backend="unknown"), "backend must"),
        (
            lambda ad, genes, group: ag.plot_embedding(
                ad,
                color=group,
                split_by=group,
                backend="matplotlib",
            ),
            "unsplit embeddings",
        ),
        (
            lambda ad, genes, group: ag.plot_embedding(ad, backend="matplotlib"),
            "point-density",
        ),
        (
            lambda ad, genes, group: ag.plot_dotplot(
                ad,
                genes,
                group,
                split_by=group,
                backend="matplotlib",
            ),
            "unsplit dotplots",
        ),
        (
            lambda ad, genes, group: ag.plot_matrixplot(
                ad,
                genes,
                group,
                split_by=group,
                backend="matplotlib",
            ),
            "unsplit matrixplots",
        ),
    ],
)
def test_fast_backend_rejects_unsupported_layouts(
    adata,
    markers,
    group_key,
    call,
    match,
):
    with pytest.raises(ValueError, match=match):
        call(adata, markers, group_key)
