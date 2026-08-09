"""Executable companion to the grammar-of-graphics vignette."""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from _fixture import make_adata
from plotnine import facet_wrap, geom_point, ggplot, scale_color_brewer, theme_classic

from ggann import aes, gganndata, obs, obsm


def main() -> None:
    adata = make_adata()
    plot = (
        gganndata(
            adata,
            aes(
                x=obsm("umap", 0),
                y=obsm("umap", 1),
                color=obs("cell_type"),
                group=obs("condition"),
            ),
            max_matrix_values=2 * adata.n_obs,
            add_theme=False,
        )
        + geom_point(size=1.8, alpha=0.85)
        + scale_color_brewer(type="qual", palette="Set2")
        + facet_wrap("condition")
        + theme_classic()
    )
    assert isinstance(plot, ggplot)
    plot.draw()
    plt.close("all")


if __name__ == "__main__":
    main()
