"""Figures for docs/index.md and docs/notes.md.

Purpose: render the home-page hero panels and the notes-page timing plot.
Run: python examples/docs_figures/index_notes.py
Inputs: data/pbmc3k_processed.h5ad (offline Scanpy dataset cache);
    timings typed from docs/notes.md (measured values, see the table there).
Outputs: docs/images/index/{umap,dotplot,violin}.png, docs/images/notes/timing.png
"""

from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")

import pandas as pd
import scanpy as sc
from plotnine import (
    aes,
    element_blank,
    element_text,
    facet_wrap,
    geom_col,
    geom_text,
    labs,
    scale_fill_manual,
    scale_y_continuous,
    theme,
)

import ggann as ag
from ggann import gganndata

ROOT = Path(__file__).resolve().parents[2]
GROUP = "louvain"
MARKERS = ["CD3D", "CD8A", "NKG7", "GNLY", "MS4A1", "FCGR3A", "CST3"]

sc.settings.datasetdir = ROOT / "data"
adata = sc.datasets.pbmc3k_processed()

# Hero panels ------------------------------------------------------------------
out = ROOT / "docs" / "images" / "index"
out.mkdir(parents=True, exist_ok=True)

umap = ag.plot_embedding(adata, "umap", color=GROUP)
dots = ag.plot_dotplot(adata, MARKERS, GROUP)
violin = ag.plot_violin(adata, MARKERS[:3], GROUP)

umap.save(out / "umap.png", width=4.5, height=4, dpi=80, verbose=False)
dots.save(out / "dotplot.png", width=5.5, height=4, dpi=80, verbose=False)
violin.save(out / "violin.png", width=6.5, height=4, dpi=80, verbose=False)

# Timing plot ------------------------------------------------------------------
# Numbers copied from docs/notes.md: matched Scanpy benchmark (20,000 x 10,000
# CSR, ms, end to end) and exact vs binned point density (s).
out = ROOT / "docs" / "images" / "notes"
out.mkdir(parents=True, exist_ok=True)

rows = [
    ("Embedding (ms)", "Scanpy", 426.7),
    ("Embedding (ms)", "ggann", 770.4),
    ("Dotplot (ms)", "Scanpy", 367.9),
    ("Dotplot (ms)", "ggann", 738.0),
    ("Matrixplot (ms)", "Scanpy", 346.5),
    ("Matrixplot (ms)", "ggann", 656.7),
    ("Point density (ms)", "exact", 27000.0),
    ("Point density (ms)", "binned", 8.0),
]
timing = pd.DataFrame(rows, columns=["panel", "method", "ms"])
timing["label"] = timing["ms"].map(lambda v: f"{v:,.0f}")
timing["method"] = pd.Categorical(timing["method"], ["Scanpy", "exact", "ggann", "binned"])

# gganndata reads an AnnData, so the table rides along as obs.
timing_ad = ad.AnnData(obs=timing.set_index(timing.index.astype(str)))
timing_plot = (
    gganndata(timing_ad, aes("method", "ms", fill="method"))
    + geom_col(width=0.7, show_legend=False)
    + geom_text(aes(label="label"), va="bottom", size=8)
    + facet_wrap("panel", scales="free", ncol=4)
    + scale_fill_manual(
        {"Scanpy": "#999999", "ggann": "#0072B2", "exact": "#999999", "binned": "#0072B2"}
    )
    + labs(y="time (ms)")
    + scale_y_continuous(expand=(0, 0, 0.12, 0))
    + theme(strip_text=element_text(size=8), axis_title_x=element_blank())
)
timing_plot.save(out / "timing.png", width=9, height=3.2, dpi=80, verbose=False)
