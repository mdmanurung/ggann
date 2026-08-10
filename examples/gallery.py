"""Generate the gallery images under ``docs/images``.

Run ``python examples/gallery.py`` with all optional plotting extras installed.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

import ggann as ag

ROOT = Path(__file__).parent.parent
OUT = ROOT / "docs" / "images"


def _save(plot, name: str, **kwargs) -> None:
    path = OUT / name
    plot.save(path, verbose=False, **kwargs)
    print("wrote", path.relative_to(Path(__file__).parent.parent))


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    sc.settings.datasetdir = ROOT / "data"
    adata = sc.datasets.pbmc3k_processed()
    # A median split of the measured library size, used by the composition and
    # trend panels in place of a treatment variable this donor sample lacks.
    adata.obs["depth"] = pd.Categorical(
        np.where(adata.obs["n_counts"] < adata.obs["n_counts"].median(), "low", "high"),
        categories=["low", "high"],
        ordered=True,
    )
    # pbmc3k_processed.X is already restricted to the highly variable genes, so
    # flagging them makes plot_correlation's default gene selection explicit.
    adata.var["highly_variable"] = True
    group = "louvain"
    markers = ["CD3D", "CD8A", "NKG7", "GNLY", "MS4A1", "FCGR3A", "CST3"]
    sc.tl.rank_genes_groups(adata, group, method="wilcoxon", n_genes=50)

    _save(
        ag.plot_embedding(adata, "umap", color=group),
        "umap_clusters.png",
        width=7,
        height=5,
        dpi=100,
    )
    _save(
        ag.plot_features(adata, markers[:4], basis="umap"),
        "features_grid.png",
        width=8,
        height=6,
        dpi=100,
    )
    _save(
        ag.plot_dotplot(adata, markers, group),
        "dotplot.png",
        width=7,
        height=5,
        dpi=100,
    )
    _save(
        ag.plot_rank_genes_dotplot(adata, n_genes=3),
        "de_dotplot.png",
        width=9,
        height=5,
        dpi=100,
    )
    _save(
        ag.plot_volcano(adata, group="NK cells"),
        "volcano.png",
        width=6.5,
        height=5,
        dpi=100,
    )
    _save(
        ag.plot_proportions(adata, group, split_by="depth", position="fill"),
        "proportions.png",
        width=6,
        height=5,
        dpi=100,
    )
    _save(
        ag.plot_stacked_violin(adata, markers, group),
        "stacked_violin.png",
        width=6,
        height=8,
        dpi=100,
    )
    _save(
        ag.plot_qc_violin(
            adata,
            metrics=["n_genes", "percent_mito", "n_counts"],
            group_by=group,
        ),
        "qc_violin.png",
        width=7,
        height=8,
        dpi=90,
    )
    _save(
        ag.plot_embedding(adata, "umap", color=group, label=True),
        "umap_labelled.png",
        width=7,
        height=5,
        dpi=100,
    )
    _save(
        ag.plot_density(adata, ["CD3D", "NKG7"], basis="umap", joint=True),
        "density.png",
        width=10,
        height=3.2,
        dpi=100,
    )
    _save(
        ag.plot_box(adata, markers[:4], group),
        "box.png",
        width=6,
        height=8,
        dpi=100,
    )
    _save(
        ag.plot_expression_bar(adata, markers[:4], group),
        "expression_bar.png",
        width=6,
        height=8,
        dpi=100,
    )
    _save(
        ag.plot_expression_line(adata, markers[:3], x="depth", group_by=group),
        "expression_line.png",
        width=7,
        height=7,
        dpi=100,
    )
    _save(
        ag.plot_correlation(adata, group, cluster=True),
        "correlation.png",
        width=6,
        height=5,
        dpi=100,
    )

    de = ag.rank_genes_df(adata)
    selected_groups = ["CD14+ Monocytes", "B cells", "NK cells", "Dendritic cells"]
    marker_sets = {
        name: list(de[de["group"] == name].head(20)["names"]) for name in selected_groups
    }
    upset_path = OUT / "upset.png"
    ag.plot_upset(marker_sets, min_cardinality=1).save(upset_path, dpi=100)
    plt.close("all")
    print("wrote", upset_path.relative_to(Path(__file__).parent.parent))


if __name__ == "__main__":
    main()
