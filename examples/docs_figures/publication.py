"""Purpose: render the figures shown on docs/publication.md.

Run: python examples/docs_figures/publication.py
Inputs:  data/pbmc3k_processed.h5ad (scanpy dataset cache)
Outputs: docs/images/publication/*.png (written by save_publication)
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import numpy as np
import pandas as pd
import scanpy as sc
from plotnine import theme

import ggann as ag

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "docs" / "images" / "publication"
OUT.mkdir(parents=True, exist_ok=True)

# Setup -----------------------------------------------------------------------
sc.settings.datasetdir = ROOT / "data"
adata = sc.datasets.pbmc3k_processed()
genes = ["CD3D", "MS4A1", "NKG7", "GNLY", "CST3"]
lineage = {
    "CD4 T cells": "T cell", "CD8 T cells": "T cell", "B cells": "B cell",
    "NK cells": "NK cell", "CD14+ Monocytes": "Monocyte",
    "FCGR3A+ Monocytes": "Monocyte", "Dendritic cells": "Dendritic",
    "Megakaryocytes": "Megakaryocyte",
}  # fmt: skip
order = ["T cell", "B cell", "NK cell", "Monocyte", "Dendritic", "Megakaryocyte"]
adata.obs["lineage"] = pd.Categorical(adata.obs["louvain"].map(lineage), order)
adata.obs["depth"] = pd.Categorical(
    np.where(adata.obs["n_counts"] < adata.obs["n_counts"].median(), "low", "high"),
    categories=["low", "high"],
)

# One figure: 89 mm wide --------------------------------------------------------
with ag.style_context("single-column"):
    umap = ag.plot_embedding(adata, "umap", color="lineage", size=1.0)
ag.save_publication(umap, OUT / "single_column.png", width="single-column", height=70, dpi=150)

# Palettes ----------------------------------------------------------------------
palette = ag.publication_palette("qualitative", categories=order)
adata.uns["lineage_colors"] = list(palette.values())
with ag.style_context("single-column"):
    dots = ag.plot_dotplot(adata, genes, group_by="lineage", use_raw=True)
ag.save_publication(dots, OUT / "palette_dotplot.png", width="single-column", height=70, dpi=150)

# Four panels, 183 mm -----------------------------------------------------------
with ag.style_context("double-column"):
    panels = [
        ag.plot_embedding(adata, "umap", color="lineage", size=1.0),
        ag.plot_dotplot(adata, genes, group_by="lineage", use_raw=True),
        ag.plot_violin(adata, ["CD3D"], group_by="lineage", use_raw=True)
        + theme(legend_position="none"),
        ag.plot_proportions(adata, group_by="lineage", split_by="depth"),
    ]
    figure = ag.compose(panels, ncol=2, widths=(0.95, 1.25), gap=2, tag_levels="auto") + theme(
        legend_position="right"
    )
ag.save_publication(figure, OUT / "four_panels.png", width="double-column", height=140, dpi=150)
