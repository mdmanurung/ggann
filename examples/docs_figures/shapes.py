"""Figures for docs/shapes.md: Bezier curves, pie, ring and radar.

Run: python examples/docs_figures/shapes.py
Outputs: docs/images/shapes/*.png
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import anndata as ad
import numpy as np
import pandas as pd
from plotnine import (
    aes,
    arrow,
    coord_equal,
    geom_point,
    geom_text,
    ggplot,
    scale_color_manual,
    theme_minimal,
)

import ggann as ag

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "docs" / "images" / "shapes"
OUT.mkdir(parents=True, exist_ok=True)

# Curves -----------------------------------------------------------------------
# One segment from (0, 0) to (4, 0) drawn at three curvatures, left to right.
CURVATURES = {"-0.5": -0.5, "0": 0.0, "0.5": 0.5}
segments = pd.DataFrame(
    {
        "x": 0.0,
        "y": 0.0,
        "xend": 4.0,
        "yend": 0.0,
        "curvature": list(CURVATURES),
    }
)
start = pd.DataFrame({"x": [0.0], "y": [0.0]})
labels = pd.DataFrame(
    {
        "x": [2.0, 2.0, 2.0],
        "y": [1.45, 0.2, -1.45],
        "label": ["curvature = -0.5", "0", "curvature = 0.5"],
    }
)
palette = {"-0.5": "#2166ac", "0": "#777777", "0.5": "#b2182b"}

plot = (
    ggplot(segments, aes("x", "y", xend="xend", yend="yend", color="curvature"))
    + ag.geom_curve(
        data=segments[segments["curvature"] == "-0.5"],
        curvature=-0.5,
        size=0.8,
        arrow=arrow(length=0.15),
    )
    + ag.geom_curve(
        data=segments[segments["curvature"] == "0"], curvature=0, size=0.8, arrow=arrow(length=0.15)
    )
    + ag.geom_curve(
        data=segments[segments["curvature"] == "0.5"],
        curvature=0.5,
        size=0.8,
        arrow=arrow(length=0.15),
    )
    + geom_point(aes("x", "y"), data=start, inherit_aes=False, size=2)
    + geom_text(aes("x", "y", label="label"), data=labels, inherit_aes=False, size=8)
    + scale_color_manual(values=palette, guide=None)
    + coord_equal()
    + theme_minimal()
)
plot.save(OUT / "curve.png", width=5, height=3, dpi=80, verbose=False)

# Pie ----------------------------------------------------------------------------
# Four cell types with weights 4, 3, 2, 1: slices run clockwise from 12 o'clock.
ORDER = ["T", "B", "NK", "Myeloid"]
CELL_TYPES = pd.DataFrame({"cell_type": ["T", "B", "NK", "Myeloid"], "n": [4, 3, 2, 1]})
pie = ag.ggpie(
    CELL_TYPES,
    "cell_type",
    "n",
    order=ORDER,
    label="percent",
    palette=["#4e79a7", "#f28e2b", "#59a14f", "#e15759"],
)
pie.save(OUT / "pie.png", width=4, height=3, dpi=80, verbose=False)

# Ring ---------------------------------------------------------------------------
# plot_proportions(kind="ring") on a table of cells: one ring per condition.
CELLS_PER_TYPE = {
    "control": {"T": 40, "B": 30, "NK": 20, "Myeloid": 10},
    "treated": {"T": 20, "B": 20, "NK": 20, "Myeloid": 40},
}
obs = pd.DataFrame(
    [
        {"cell_type": cell_type, "condition": condition}
        for condition, counts in CELLS_PER_TYPE.items()
        for cell_type, n in counts.items()
        for _ in range(n)
    ]
)
obs["cell_type"] = pd.Categorical(obs["cell_type"], categories=ORDER)
obs["condition"] = pd.Categorical(obs["condition"], categories=list(CELLS_PER_TYPE))
obs.index = [f"cell_{i}" for i in range(len(obs))]
adata = ad.AnnData(np.zeros((len(obs), 1)), obs=obs)
adata.uns["cell_type_colors"] = ["#4e79a7", "#f28e2b", "#59a14f", "#e15759"]
ring = ag.plot_proportions(adata, "cell_type", "condition", kind="ring")
ring.save(OUT / "ring.png", width=5, height=3, dpi=80, verbose=False)

# Radar --------------------------------------------------------------------------
# Five markers scored in two groups on one radius axis that starts at 0.
SCORES = pd.DataFrame(
    {
        "marker": ["CD3", "CD19", "CD56", "CD14", "CD34"] * 2,
        "score": [4, 2, 3, 1, 2, 1, 3, 2, 4, 3],
        "group": ["control"] * 5 + ["treated"] * 5,
    }
)
radar = ag.ggradar(
    SCORES,
    "marker",
    "score",
    series="group",
    axis_order=["CD3", "CD19", "CD56", "CD14", "CD34"],
    palette={"control": "#4e79a7", "treated": "#e15759"},
)
radar.save(OUT / "radar.png", width=4.5, height=3.5, dpi=80, verbose=False)
