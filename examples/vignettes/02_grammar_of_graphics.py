"""Test a condition-associated expression pattern with the plotnine grammar."""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from _fixture import make_adata
from plotnine import facet_wrap, geom_point, ggplot, labs, scale_color_cmap, theme_classic

from ggann import aes, gene, gganndata, obs, obsm


def main() -> None:
    adata = make_adata()
    plot = (
        gganndata(
            adata,
            aes(
                x=obsm("umap", 0),
                y=obsm("umap", 1),
                color=gene("MKI67", layer="logcounts"),
                group=obs("condition"),
            ),
            max_matrix_values=3 * adata.n_obs,
            add_theme=False,
        )
        + geom_point(size=1.8, alpha=0.85)
        + scale_color_cmap(cmap_name="viridis")
        + facet_wrap("condition")
        + labs(title="MKI67 across conditions", color="log expression")
        + theme_classic()
    )
    assert isinstance(plot, ggplot)
    condition_means = plot.data.groupby("condition", observed=True)["MKI67"].mean()
    assert condition_means["stimulated"] > condition_means["control"]
    plot.draw()
    plt.close("all")


if __name__ == "__main__":
    main()
