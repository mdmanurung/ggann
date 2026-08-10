"""Executable comparison of ggann + annplyr and standard Scanpy plotting."""

from __future__ import annotations

import json
import tempfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from _fixture import CELL_TYPES, MARKERS, fingerprint, load_adata
from matplotlib.figure import Figure
from plotnine import aes, facet_wrap, geom_point, scale_color_cmap, theme_classic

import ggann as ag
from ggann._aggregate import aggregate_expression

FIGURE_SIZE = (6.0, 4.5)
DPI = 80
GENES = MARKERS
ROOT = Path(__file__).resolve().parents[2]
BENCHMARK = ROOT / "benchmarks" / "results" / "scanpy-extended-csr.json"


def _scanpy_figure(output) -> Figure:
    if isinstance(output, Figure):
        return output
    if hasattr(output, "make_figure") and getattr(output, "fig", None) is None:
        output.make_figure()
    figure = getattr(output, "fig", None)
    if figure is None:
        raise TypeError(f"Could not recover a Scanpy figure from {type(output)!r}")
    return figure


def _save_ggann(plot, path: Path) -> None:
    figure = plot.draw(show=False)
    figure.set_size_inches(*FIGURE_SIZE, forward=True)
    figure.savefig(path, format="png", dpi=DPI)
    plt.close(figure)


def _save_scanpy(output, path: Path) -> None:
    figure = _scanpy_figure(output)
    figure.set_size_inches(*FIGURE_SIZE, forward=True)
    figure.savefig(path, format="png", dpi=DPI)
    plt.close(figure)


def _validate_prepared_payloads(adata) -> None:
    grammar = ag.gganndata(
        adata,
        aes(
            x=ag.obsm("umap", 0),
            y=ag.obsm("umap", 1),
            color=ag.gene("NKG7", layer="logcounts"),
        ),
        max_matrix_values=3 * adata.n_obs,
        add_theme=False,
    )
    expected = sc.get.obs_df(adata, keys=["NKG7"], layer="logcounts", use_raw=False)
    np.testing.assert_allclose(grammar.data["UMAP_1"], adata.obsm["X_umap"][:, 0])
    np.testing.assert_allclose(grammar.data["UMAP_2"], adata.obsm["X_umap"][:, 1])
    np.testing.assert_allclose(grammar.data["NKG7"], expected["NKG7"])

    ggann_summary = aggregate_expression(
        adata,
        GENES,
        "louvain",
        use_raw=False,
    )
    wide = sc.get.obs_df(adata, keys=["louvain", *GENES], use_raw=False).set_index("louvain")
    scanpy_means = wide.groupby(level=0, observed=True, sort=False).mean().reindex(CELL_TYPES)
    scanpy_fraction = (
        (wide > 0).groupby(level=0, observed=True, sort=False).mean().reindex(CELL_TYPES)
    )
    ggann_means = ggann_summary.pivot(
        index="louvain", columns="feature", values="mean_expression"
    ).reindex(index=CELL_TYPES, columns=GENES)
    ggann_fraction = ggann_summary.pivot(
        index="louvain", columns="feature", values="fraction"
    ).reindex(index=CELL_TYPES, columns=GENES)
    np.testing.assert_allclose(ggann_means, scanpy_means[GENES], rtol=1e-6, atol=1e-7)
    np.testing.assert_allclose(ggann_fraction, scanpy_fraction[GENES], rtol=1e-6, atol=1e-7)


def main() -> None:
    adata = load_adata(storage="csr")
    before = fingerprint(adata)
    _validate_prepared_payloads(adata)

    benchmark = json.loads(BENCHMARK.read_text())
    assert benchmark["benchmark_kind"] == "ggann_scanpy_matched"
    assert benchmark["metadata"]["repeats"] >= 5
    assert benchmark["release_gates"]["status"] == "fail"

    ggann_plots = {
        "embedding-categorical": ag.plot_embedding(
            adata,
            "umap",
            color="louvain",
            pointdensity=False,
        ),
        "embedding-gene": ag.plot_embedding(
            adata,
            "umap",
            color=ag.gene("NKG7", layer="logcounts"),
            pointdensity=False,
        ),
        "dotplot": ag.plot_dotplot(
            adata,
            GENES,
            group_by="louvain",
            use_raw=False,
            categories_order=CELL_TYPES,
        ),
        "matrixplot": ag.plot_matrixplot(
            adata,
            GENES,
            group_by="louvain",
            use_raw=False,
            categories_order=CELL_TYPES,
        ),
    }
    scanpy_plots = {
        "embedding-categorical": sc.pl.embedding(
            adata,
            "umap",
            color="louvain",
            show=False,
            return_fig=True,
        ),
        "embedding-gene": sc.pl.embedding(
            adata,
            "umap",
            color="NKG7",
            layer="logcounts",
            use_raw=False,
            show=False,
            return_fig=True,
        ),
        "dotplot": sc.pl.dotplot(
            adata,
            GENES,
            "louvain",
            use_raw=False,
            categories_order=CELL_TYPES,
            figsize=FIGURE_SIZE,
            show=False,
            return_fig=True,
        ),
        "matrixplot": sc.pl.matrixplot(
            adata,
            GENES,
            "louvain",
            use_raw=False,
            categories_order=CELL_TYPES,
            figsize=FIGURE_SIZE,
            show=False,
            return_fig=True,
        ),
    }

    with tempfile.TemporaryDirectory(prefix="ggann-scanpy-vignette-") as temporary:
        output = Path(temporary)
        for name, plot in ggann_plots.items():
            _save_ggann(plot, output / f"ggann-{name}.png")
        for name, plot in scanpy_plots.items():
            _save_scanpy(plot, output / f"scanpy-{name}.png")
        images = sorted(output.glob("*.png"))
        assert len(images) == 8
        assert all(path.stat().st_size > 1_000 for path in images)

        publication_spec = (
            ag.gganndata(
                adata,
                aes(
                    x=ag.obsm("umap", 0),
                    y=ag.obsm("umap", 1),
                    color=ag.gene("NKG7", layer="logcounts"),
                    group=ag.obs("compartment"),
                ),
                max_matrix_values=3 * adata.n_obs,
                add_theme=False,
            )
            + geom_point(size=1.8, alpha=0.85)
            + facet_wrap("compartment")
            + scale_color_cmap(cmap_name="magma")
            + theme_classic()
        )
        _save_ggann(publication_spec, output / "ggann-publication-spec.png")

    assert isinstance(ggann_plots["dotplot"].data, pd.DataFrame)
    assert fingerprint(adata) == before
    plt.close("all")


if __name__ == "__main__":
    main()
