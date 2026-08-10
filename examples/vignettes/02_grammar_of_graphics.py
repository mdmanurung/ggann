"""Test a lineage-associated expression pattern with the plotnine grammar."""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from _fixture import load_adata
from plotnine import facet_wrap, geom_point, ggplot, labs, scale_color_cmap, theme_classic

from ggann import aes, gene, gganndata, obs, obsm


def main() -> None:
    adata = load_adata()
    plot = (
        gganndata(
            adata,
            aes(
                x=obsm("umap", 0),
                y=obsm("umap", 1),
                color=gene("NKG7", layer="logcounts"),
                group=obs("compartment"),
            ),
            max_matrix_values=3 * adata.n_obs,
            add_theme=False,
        )
        + geom_point(size=0.8, alpha=0.85)
        + scale_color_cmap(cmap_name="viridis")
        + facet_wrap("compartment")
        + labs(title="NKG7 across PBMC lineage compartments", color="log expression")
        + theme_classic()
    )
    assert isinstance(plot, ggplot)
    compartment_means = plot.data.groupby("compartment", observed=True)["NKG7"].mean()
    assert compartment_means["Lymphoid"] > compartment_means["Myeloid"]
    plot.draw()
    plt.close("all")


if __name__ == "__main__":
    main()
