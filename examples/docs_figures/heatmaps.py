"""Figures for docs/heatmaps.md: clustered heatmaps and UpSet (grid backends).

Run: python examples/docs_figures/heatmaps.py
Outputs: docs/images/heatmaps/*.png
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

import ggann as ag

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "docs" / "images" / "heatmaps"
OUT.mkdir(parents=True, exist_ok=True)

GROUP = "louvain"
MARKERS = ["CD3D", "CD8A", "NKG7", "GNLY", "MS4A1", "FCGR3A", "CST3"]

sc.settings.datasetdir = ROOT / "data"
adata = sc.datasets.pbmc3k_processed()
adata.obs["depth"] = pd.Categorical(
    np.where(adata.obs["n_counts"] < adata.obs["n_counts"].median(), "low", "high"),
    categories=["low", "high"],
    ordered=True,
)

short = {
    "CD4 T cells": "CD4 T",
    "CD8 T cells": "CD8 T",
    "CD14+ Monocytes": "CD14 Mono",
    "FCGR3A+ Monocytes": "FCGR3A Mono",
    "Dendritic cells": "DC",
    "Megakaryocytes": "Mk",
    "NK cells": "NK",
    "B cells": "B",
}
# Short labels in a slim copy that has no other per-group obs column
short_ad = adata.copy()
short_ad.obs = pd.DataFrame({"cell_type": adata.obs[GROUP].map(short).astype("category")})

# Grouped means --------------------------------------------------------------
plt.figure(figsize=(5, 3.5))
cm = ag.plot_clustermap(
    adata,
    MARKERS,
    group_by=GROUP,
    standard_scale="var",
    row_dendrogram=True,
    col_dendrogram=True,
)
plt.savefig(OUT / "clustermap_groups.png", dpi=80, bbox_inches="tight")
plt.close("all")

# Single cells with annotation bars ------------------------------------------
cells = adata[adata.obs.sample(300, random_state=0).index]
plt.figure(figsize=(6, 3.5))
cm = ag.plot_clustermap(
    cells,
    MARKERS,
    annotations=[GROUP, "depth"],
    z_score=0,
    show_colnames=False,
)
plt.savefig(OUT / "clustermap_cells.png", dpi=80, bbox_inches="tight")
plt.close("all")

# tidy_heatmap from AnnData --------------------------------------------------
heatmap = (
    ag.tidy_heatmap(
        short_ad,
        genes=MARKERS,
        group_by="cell_type",
        scale="row",
        show_colnames=False,  # the tile strip names the groups
    )
    .annotation_tile("cell_type")
    .layer_point("mean_expression > 1", size=2)
)
plt.close("all")
heatmap.save(OUT / "tidy_adata", width="double-column", height=90, formats=("png",), dpi=80)
plt.close("all")

# tidy_heatmap from a long table ---------------------------------------------
rng = np.random.default_rng(0)
lineage = {"T": "lymphoid", "B": "lymphoid", "NK": "lymphoid", "Mono": "myeloid", "DC": "myeloid"}
df = pd.DataFrame(
    [
        {
            "feature": f"G{i:02d}",
            "cell_type": cell_type,
            "mean": rng.gamma(2.0, 1.0),
            "padj": rng.random() ** 2,
            "lineage": lineage[cell_type],
            "n_cells": 20 * (k + 1),
        }
        for i in range(14)
        for k, cell_type in enumerate(lineage)
    ]
)
long_heatmap = (
    ag.tidy_heatmap(df, row="feature", column="cell_type", value="mean", scale="row")
    .annotation_tile("lineage")
    .annotation_bar("n_cells")
    .split_rows(3)
    .layer_point("padj < 0.05")
)
plt.close("all")
long_heatmap.save(OUT / "tidy_long", width=150, height=100, formats=("png",), dpi=80)
plt.close("all")

# Export: auto-fit to a single column ---------------------------------------
plt.figure(figsize=(90 / 25.4, 70 / 25.4))
cm = ag.plot_clustermap(short_ad, MARKERS, group_by="cell_type", standard_scale="var", legend=False)
ag.save_publication(
    cm, OUT / "export_single_column", width="single-column", height=70, formats=("png",), dpi=150
)
plt.close("all")

# UpSet: top DE genes per cluster -------------------------------------------
sc.tl.rank_genes_groups(adata, GROUP, method="wilcoxon", n_genes=50, key_added="rank_wilcoxon")
de = ag.rank_genes_df(adata, key="rank_wilcoxon")
groups = ["CD4 T cells", "CD8 T cells", "NK cells", "B cells"]
sets = {g: de[de["group"] == g].head(40)["names"] for g in groups}
plt.figure(figsize=(6, 3.5))
up = ag.plot_upset(sets)
plt.savefig(OUT / "upset.png", dpi=80, bbox_inches="tight")
plt.close("all")
