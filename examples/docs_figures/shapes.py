"""Figures for docs/shapes.md: Bezier curves.

Run: python examples/docs_figures/shapes.py
Outputs: docs/images/shapes/*.png
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

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
