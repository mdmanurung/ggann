"""Figures for docs/grammar.md. Each block is the code shown on that page.

Run: python examples/docs_figures/grammar.py
Inputs:  data/pbmc3k_processed.h5ad
Outputs: docs/images/grammar/*.png
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import annplyr as ap
import numpy as np
import scanpy as sc
from annplyr.data import obs, raw
from plotnine import (
    aes,
    element_text,
    geom_boxplot,
    geom_col,
    geom_jitter,
    geom_point,
    geom_violin,
    ggplot,
    scale_y_continuous,
    theme,
)

import ggann as ag
from ggann import gene, gganndata, obsm

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "docs" / "images" / "grammar"
OUT.mkdir(parents=True, exist_ok=True)


def save(plot, name, width=5.5, height=4):
    plot.save(OUT / f"{name}.png", width=width, height=height, dpi=80, verbose=False)


sc.settings.datasetdir = ROOT / "data"
adata = sc.datasets.pbmc3k_processed()

# --- 1. gganndata in one example ------------------------------------------- #
plot = gganndata(adata, aes("UMAP_1", "UMAP_2", color="louvain")) + geom_point()
save(plot, "first_scatter")

# --- 2. what names resolve ------------------------------------------------- #
plot = (
    gganndata(
        adata,
        aes(
            x=obsm("umap", 0),
            y=obsm("umap", 1),
            color=gene("CD3D", use_raw=True),
        ),
    )
    + geom_point()
    + ag.theme_ggann()
)
save(plot, "resolve", width=6.5)

# --- 3. helper vs grammar -------------------------------------------------- #
helper = ag.plot_violin(adata, ["NKG7"], "louvain", use_raw=True, add_box=False)
save(helper, "violin_helper")

grammar = (
    gganndata(adata, aes("louvain", gene("NKG7", use_raw=True)))
    + geom_violin(scale="width")
    + ag.theme_ggann()
    + theme(axis_text_x=element_text(rotation=45, ha="right"))
)
save(grammar, "violin_grammar")

extended = helper + ag.stat_central_tendency()
save(extended, "violin_extended")

# --- 4. statistics layers -------------------------------------------------- #
pair = adata[adata.obs["louvain"].isin(["CD4 T cells", "CD14+ Monocytes"])].copy()
pair.obs["louvain"] = pair.obs["louvain"].cat.remove_unused_categories()
stats_plot = (
    gganndata(pair, aes("louvain", gene("CD3D", use_raw=True)))
    + geom_violin()
    + ag.stat_compare_means(method="wilcox.test")
    + ag.theme_ggann()
)
save(stats_plot, "stat_compare", width=4)

pwc = (
    gganndata(
        adata[adata.obs["louvain"].isin(["CD4 T cells", "CD8 T cells", "NK cells"])],
        aes("louvain", gene("GNLY", use_raw=True)),
    )
    + geom_boxplot(outlier_shape="")
    + ag.stat_pwc(method="wilcox.test", label="p.adj.signif")
    + scale_y_continuous(expand=(0.05, 0, 0.2, 0))
    + ag.theme_ggann()
)
save(pwc, "stat_pwc", width=4)

cor = (
    gganndata(adata, aes(gene("CD3D", use_raw=True), gene("CD8A", use_raw=True)))
    + geom_jitter(alpha=0.3, size=1, width=0.05, height=0.05)
    + ag.stat_cor(method="spearman", label_x_npc=0.4)
    + ag.theme_ggann()
)
save(cor, "stat_cor", width=4)

# --- 5. pseudobulk --------------------------------------------------------- #
# pbmc3k is one donor. Four random pseudo-donors stand in, so the plot shows the
# layout, not biology. .raw is log-normalised, so use mean mode and skip_checks.
rng = np.random.default_rng(0)
adata.obs["donor"] = rng.choice(["D1", "D2", "D3", "D4"], adata.n_obs)
pb = ag.pseudobulk(
    adata,
    sample_col="donor",
    group_by="louvain",
    use_raw=True,
    mode="mean",
    min_cells=5,
    skip_checks=True,
)
plot = (
    gganndata(pb, aes("louvain", gene("NKG7", use_raw=False), color="donor"))
    + geom_jitter(width=0.15, height=0, size=2)
    + ag.theme_ggann()
    + theme(axis_text_x=element_text(rotation=45, ha="right"))
)
save(plot, "pseudobulk")

# --- 6. annplyr pipelines -------------------------------------------------- #
tidy = ap.tidy(adata)
scored = tidy.filter(obs.percent_mito < 0.03).mutate(
    t_score=(raw.CD3D + raw.CD3E) / 2,
)
plot = gganndata(scored, aes("UMAP_1", "UMAP_2", color="t_score")) + geom_point() + ag.theme_ggann()
save(plot, "annplyr_score")

summary = tidy.summarize(
    mean_NKG7=ap.mean(raw.NKG7),
    by=obs.louvain,
    max_matrix_values=adata.n_obs,
)
plot = (
    ggplot(summary, aes("louvain", "mean_NKG7"))
    + geom_col()
    + ag.theme_ggann()
    + theme(axis_text_x=element_text(rotation=45, ha="right"))
)
save(plot, "annplyr_summary")

# --- 7. large data --------------------------------------------------------- #
plot = (
    gganndata(
        adata,
        aes("UMAP_1", "UMAP_2", color=gene("CD3D", use_raw=True)),
        max_matrix_values=3 * adata.n_obs,
    )
    + geom_point()
)
