"""Render scanpy and ggann plots side by side.

For each ``sc.pl.*`` figure we render scanpy's output and the ggann equivalent to
``docs/images/scanpy/<name>_scanpy.png`` / ``<name>_ggann.png`` on the same
``pbmc3k_processed`` data for the comparisons in ``docs/scanpy_parity.md``.

Run: ``python examples/reproduce_scanpy.py``.
"""

from __future__ import annotations

import os

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

import ggann as ag

GROUP = "louvain"
MARKERS = ["CD3D", "CD8A", "NKG7", "GNLY", "MS4A1", "FCGR3A", "CST3"]


def _adata():
    sc.settings.datasetdir = os.path.join(os.path.dirname(__file__), "..", "data")
    adata = sc.datasets.pbmc3k_processed()
    # A median split of the measured library size, standing in for a treatment
    # variable this single-donor sample does not have.
    adata.obs["depth"] = pd.Categorical(
        np.where(adata.obs["n_counts"] < adata.obs["n_counts"].median(), "low", "high"),
        categories=["low", "high"],
        ordered=True,
    )
    # pbmc3k_processed.X is already restricted to the highly variable genes.
    adata.var["highly_variable"] = True
    sc.tl.rank_genes_groups(adata, GROUP, method="wilcoxon", n_genes=50)
    sc.tl.dendrogram(adata, groupby=GROUP)
    return adata


def _save_scanpy(fn, path):
    plt.close("all")
    fn()
    plt.savefig(path, dpi=80, bbox_inches="tight")
    plt.close("all")


def _save_ggann(plot, path):
    plot.save(path, width=5.5, height=4, dpi=80, verbose=False)


def pairs(adata):
    return {
        "umap": (
            lambda: sc.pl.umap(adata, color=GROUP, show=False),
            lambda: ag.plot_embedding(adata, "umap", color=GROUP, label=True),
        ),
        "umap_gene": (
            lambda: sc.pl.umap(adata, color=MARKERS[:4], show=False),
            lambda: ag.plot_features(adata, MARKERS[:4], basis="umap"),
        ),
        "dotplot": (
            lambda: sc.pl.dotplot(adata, MARKERS, groupby=GROUP, show=False),
            lambda: ag.plot_dotplot(adata, MARKERS, GROUP),
        ),
        "matrixplot": (
            lambda: sc.pl.matrixplot(adata, MARKERS, groupby=GROUP, show=False),
            lambda: ag.plot_matrixplot(adata, MARKERS, GROUP, standard_scale="var"),
        ),
        "stacked_violin": (
            lambda: sc.pl.stacked_violin(adata, MARKERS, groupby=GROUP, show=False),
            lambda: ag.plot_stacked_violin(adata, MARKERS, GROUP),
        ),
        "violin": (
            lambda: sc.pl.violin(adata, MARKERS[:3], groupby=GROUP, show=False),
            lambda: ag.plot_violin(adata, MARKERS[:3], GROUP),
        ),
        "tracksplot": (
            lambda: sc.pl.tracksplot(adata, MARKERS, groupby=GROUP, show=False),
            lambda: ag.plot_tracksplot(adata, MARKERS, GROUP),
        ),
        "correlation": (
            lambda: sc.pl.correlation_matrix(adata, groupby=GROUP, show=False),
            lambda: ag.plot_correlation(adata, GROUP),
        ),
        "highest_expr_genes": (
            lambda: sc.pl.highest_expr_genes(adata, n_top=15, show=False),
            lambda: ag.plot_highest_expr_genes(adata, n=15, use_raw=True),
        ),
        "rank_genes_dotplot": (
            lambda: sc.pl.rank_genes_groups_dotplot(adata, n_genes=3, show=False),
            lambda: ag.plot_rank_genes_dotplot(adata, n_genes=3),
        ),
    }


def main():
    adata = _adata()
    out = os.path.join(os.path.dirname(__file__), "..", "docs", "images", "scanpy")
    os.makedirs(out, exist_ok=True)
    failures = []
    for name, (scfn, ggfn) in pairs(adata).items():
        try:
            _save_scanpy(scfn, os.path.join(out, f"{name}_scanpy.png"))
            _save_ggann(ggfn(), os.path.join(out, f"{name}_ggann.png"))
            print("wrote", name)
        except Exception as error:
            failures.append((name, error))
            print("FAILED", name, type(error).__name__, str(error)[:80])

    if failures:
        names = ", ".join(name for name, _ in failures)
        raise RuntimeError(f"Failed to render: {names}")


if __name__ == "__main__":
    main()
