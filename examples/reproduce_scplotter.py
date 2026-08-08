"""Render ggann counterparts to scplotter vignette figures.

scplotter's vignettes run on a Seurat ``pancreas_sub`` object; here we reproduce
the chart types and options on ``pbmc68k_reduced``. Capabilities that require
RNA velocity, lineage fits, neighbour graphs, or three-dimensional rendering
are outside this example.

Run: ``python examples/reproduce_scplotter.py`` -> writes docs/images/scplotter/.
"""

from __future__ import annotations

import os

import matplotlib

matplotlib.use("Agg")

import numpy as np
import scanpy as sc
from plotnine import (
    aes,
    coord_flip,
    element_blank,
    geom_bin2d,
    geom_density_2d,
    ggplot,
    labs,
    scale_fill_cmap,
    stat_summary,
    theme,
)

import ggann as ag

GROUP = "bulk_labels"
MARKERS = ["CD3D", "CD8A", "NKG7", "GNLY", "MS4A1", "FCGR3A", "CST3"]


def _adata():
    adata = sc.datasets.pbmc68k_reduced()
    sc.tl.rank_genes_groups(adata, GROUP, method="wilcoxon", n_genes=50)
    return adata


# name -> (builder, "scplotter call it reproduces")
def build(adata):
    coords = ag.embedding_coords(adata, "umap")
    x, y = coords.columns[:2]
    figs = {}

    # --- CellDimPlot -------------------------------------------------------- #
    figs["celldim_basic"] = (
        ag.plot_embedding(adata, "umap", color=GROUP),
        'CellDimPlot(group_by=..., reduction="UMAP")',
    )
    figs["celldim_label"] = (
        ag.plot_embedding(adata, "umap", color=GROUP, label=True),
        "CellDimPlot(..., label = TRUE)",
    )
    figs["celldim_split"] = (
        ag.plot_embedding(adata, "umap", color=GROUP, split_by="phase"),
        'CellDimPlot(..., split_by = "Phase")',
    )
    figs["celldim_feature"] = (
        ag.plot_embedding(adata, "umap", color="CD3D"),
        'FeatureStatPlot(plot_type = "dim", features = ...)',
    )
    figs["celldim_features_grid"] = (
        ag.plot_features(adata, MARKERS[:4], basis="umap"),
        'FeatureStatPlot(plot_type = "dim", features = c(...))',
    )
    # add_density = TRUE  -> compose a density-contour layer
    figs["celldim_density"] = (
        ag.plot_embedding(adata, "umap", color=GROUP)
        + geom_density_2d(
            aes(x, y), inherit_aes=False, color="black", size=0.25, alpha=0.5
        ),
        "CellDimPlot(..., add_density = TRUE)  [composed: + geom_density_2d]",
    )
    # hex = TRUE  -> compose binned counts (geom_hex-equivalent)
    figs["celldim_hex"] = (
        ggplot(coords, aes(x, y))
        + geom_bin2d(bins=28)
        + scale_fill_cmap(cmap_name="magma")
        + labs(fill="count")
        + ag.theme_ggann()
        + theme(axis_text=element_blank(), axis_ticks=element_blank()),
        "CellDimPlot(..., hex = TRUE)  [composed: + geom_bin2d]",
    )

    # --- FeatureStatPlot ---------------------------------------------------- #
    figs["feat_violin"] = (
        ag.plot_violin(adata, ["CD3D", "NKG7"], GROUP),
        "FeatureStatPlot(features = c(...), ident = ...)  # violin",
    )
    figs["feat_violin_points"] = (
        ag.plot_violin(adata, ["CD3D"], GROUP, add_points=True),
        "FeatureStatPlot(..., add_point = TRUE)",
    )
    figs["feat_box"] = (
        ag.plot_box(adata, ["CD3D", "NKG7"], GROUP),
        'FeatureStatPlot(..., plot_type = "box")',
    )
    figs["feat_bar"] = (
        ag.plot_expression_bar(adata, ["CD3D", "NKG7"], GROUP),
        'FeatureStatPlot(..., plot_type = "bar")',
    )
    figs["feat_stacked_violin"] = (
        ag.plot_stacked_violin(adata, MARKERS, GROUP),
        "FeatureStatPlot(features = c(...), stack = TRUE)  # stacked violin",
    )
    figs["feat_stats"] = (
        ag.plot_violin(adata, ["CD3D"], GROUP, stats=True),
        "FeatureStatPlot(..., comparisons = TRUE)",
    )
    # add_stat = mean -> compose a mean crossbar
    figs["feat_addstat"] = (
        ag.plot_violin(adata, ["CD3D"], GROUP, add_box=False)
        + stat_summary(fun_y=np.mean, geom="point", size=3, color="black"),
        "FeatureStatPlot(..., add_stat = mean)  [composed: + stat_summary]",
    )
    figs["feat_heatmap"] = (
        ag.plot_matrixplot(adata, MARKERS, GROUP, standard_scale="var"),
        'FeatureStatPlot(..., plot_type = "heatmap")',
    )
    figs["feat_dot"] = (
        ag.plot_dotplot(adata, MARKERS, GROUP),
        'FeatureStatPlot(..., plot_type = "dot")',
    )
    figs["feat_cor"] = (
        ag.plot_correlation(adata, GROUP, genes=MARKERS, cluster=False, annotate=True),
        'FeatureStatPlot(..., plot_type = "cor")',
    )
    # box flipped -> compose coord_flip
    figs["feat_box_flip"] = (
        ag.plot_box(adata, ["CD3D"], GROUP) + coord_flip(),
        'FeatureStatPlot(..., plot_type = "box", flip = TRUE)  [composed: + coord_flip]',
    )
    return figs


def main():
    adata = _adata()
    out = os.path.join(os.path.dirname(__file__), "..", "docs", "images", "scplotter")
    os.makedirs(out, exist_ok=True)
    for name, (plot, _call) in build(adata).items():
        plot.save(
            os.path.join(out, f"{name}.png"), width=5.5, height=4, dpi=80, verbose=False
        )
        print("wrote", name)


if __name__ == "__main__":
    main()
