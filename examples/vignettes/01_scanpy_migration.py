"""Executable companion to the Scanpy migration vignette."""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import scanpy as sc
from _fixture import make_adata
from plotnine import ggplot

import ggann as ag


def main() -> None:
    adata = make_adata()

    # Scanpy and ggann consume the same immutable AnnData object.
    sc.pl.embedding(adata, basis="umap", color="cell_type", show=False)
    plt.close("all")

    plot = ag.plot_embedding(adata, basis="umap", color="cell_type", label=True)
    assert isinstance(plot, ggplot)
    plot.draw()
    plt.close("all")


if __name__ == "__main__":
    main()
