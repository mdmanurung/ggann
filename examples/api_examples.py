"""Render one figure per public plotting function for the API docs.

Writes ``docs/images/api/ggann.<func>.png`` for each helper; the Sphinx extension
``docs/extensions/api_examples.py`` then injects the matching image into that
function's API-reference page. Run locally with all optional dependencies
installed and commit the PNGs. The docs build does not execute plotting code.

Run: ``python examples/api_examples.py``.
"""

from __future__ import annotations

import os

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import plotnine as p9
import scanpy as sc
from plotnine import geom_point

import ggann as ag
from ggann import aes, gganndata

GROUP = "louvain"
MARKERS = ["CD3D", "CD8A", "NKG7", "GNLY", "MS4A1", "FCGR3A", "CST3"]
GENE_GROUPS = {
    "T": ["CD3D", "CD8A"],
    "NK": ["NKG7", "GNLY"],
    "B": ["MS4A1"],
    "Mono": ["CST3"],
}


def _adata():
    sc.settings.datasetdir = os.path.join(os.path.dirname(__file__), "..", "data")
    adata = sc.datasets.pbmc3k_processed()
    # A median split of the measured library size, used by the composition and
    # trend panels in place of a treatment variable this donor sample lacks.
    import numpy as np
    import pandas as pd

    adata.obs["depth"] = pd.Categorical(
        np.where(adata.obs["n_counts"] < adata.obs["n_counts"].median(), "low", "high"),
        categories=["low", "high"],
        ordered=True,
    )
    # pbmc3k_processed.X is already restricted to the highly variable genes, so
    # flagging them makes plot_correlation's default gene selection explicit.
    adata.var["highly_variable"] = True
    sc.tl.rank_genes_groups(adata, GROUP, method="wilcoxon", n_genes=50)
    return adata


def _ma_frame(adata, group="NK cells"):
    """A PyDESeq2-shaped DE table for the MA-plot example, built from the marker
    result plus each gene's mean log-normalized expression (the "base mean")."""
    import numpy as np

    raw = adata.raw.to_adata() if adata.raw is not None else adata
    means = np.asarray(raw.X.mean(axis=0)).ravel()
    mean_by_gene = dict(zip(raw.var_names, means))
    df = ag.rank_genes_df(adata, group=group).rename(
        columns={"logfoldchanges": "log2FoldChange", "pvals_adj": "padj"}
    )
    df["baseMean"] = df["names"].map(mean_by_gene)
    return df.dropna(subset=["baseMean"]).set_index("names")


def _examples(adata):
    de = ag.rank_genes_df(adata)
    sel = ["CD14+ Monocytes", "B cells", "NK cells", "Dendritic cells"]
    marker_sets = {g: list(de[de["group"] == g].head(20)["names"]) for g in sel}
    return {
        "ggann.gganndata": lambda: (
            gganndata(adata, aes("UMAP_1", "UMAP_2", color=GROUP))
            + geom_point(size=1.2)
            + ag.theme_ggann()
        ),
        "ggann.plot_embedding": lambda: ag.plot_embedding(adata, "umap", color=GROUP, label=True),
        "ggann.plot_features": lambda: ag.plot_features(adata, MARKERS[:4], basis="umap"),
        "ggann.plot_density": lambda: ag.plot_density(adata, ["CD3D", "NKG7"], joint=True),
        "ggann.plot_embedding_density": lambda: ag.plot_embedding_density(adata, "umap", GROUP),
        "ggann.plot_dotplot": lambda: ag.plot_dotplot(adata, MARKERS, GROUP),
        "ggann.plot_dotplot_grouped": lambda: ag.plot_dotplot_grouped(adata, GENE_GROUPS, GROUP),
        "ggann.plot_matrixplot": lambda: ag.plot_matrixplot(
            adata, MARKERS, GROUP, standard_scale="var"
        ),
        "ggann.plot_matrixplot_grouped": lambda: ag.plot_matrixplot_grouped(
            adata, GENE_GROUPS, GROUP
        ),
        "ggann.plot_heatmap": lambda: ag.plot_heatmap(
            adata, MARKERS, GROUP, use_raw=True, standard_scale="var"
        ),
        "ggann.plot_violin": lambda: ag.plot_violin(adata, MARKERS[:3], GROUP),
        "ggann.plot_ridge": lambda: ag.plot_ridge(adata, MARKERS[:3], GROUP),
        "ggann.plot_stacked_violin": lambda: ag.plot_stacked_violin(adata, MARKERS, GROUP),
        "ggann.plot_tracksplot": lambda: ag.plot_tracksplot(adata, MARKERS, GROUP),
        "ggann.plot_dendrogram": lambda: ag.plot_dendrogram(adata, GROUP),
        "ggann.plot_box": lambda: ag.plot_box(adata, MARKERS[:3], GROUP),
        "ggann.plot_sina": lambda: ag.plot_sina(adata, MARKERS[:3], GROUP, use_raw=True),
        "ggann.plot_expression_bar": lambda: ag.plot_expression_bar(adata, MARKERS[:3], GROUP),
        "ggann.plot_expression_line": lambda: ag.plot_expression_line(
            adata, ["CD3D"], x="depth", group_by=GROUP
        ),
        "ggann.plot_proportions": lambda: ag.plot_proportions(adata, GROUP, split_by="depth"),
        "ggann.plot_correlation": lambda: ag.plot_correlation(adata, GROUP),
        "ggann.plot_rank_genes_dotplot": lambda: ag.plot_rank_genes_dotplot(adata, n_genes=3),
        "ggann.plot_rank_genes_matrixplot": lambda: ag.plot_rank_genes_matrixplot(adata, n_genes=3),
        "ggann.plot_volcano": lambda: ag.plot_volcano(adata, group="NK cells"),
        "ggann.plot_ma": lambda: ag.plot_ma(_ma_frame(adata), label_top=8),
        "ggann.plot_qc_violin": lambda: ag.plot_qc_violin(
            adata, metrics=["n_genes", "percent_mito", "n_counts"], group_by=GROUP
        ),
        "ggann.plot_qc_scatter": lambda: ag.plot_qc_scatter(
            adata, x="n_counts", y="n_genes", color=GROUP
        ),
        # use_raw: pbmc3k_processed.X is scaled (z-scored); .raw is log-normalized,
        # so "% of total counts" is meaningful there rather than blowing up.
        "ggann.plot_highest_expr_genes": lambda: ag.plot_highest_expr_genes(
            adata, n=20, use_raw=True
        ),
        "ggann.plot_variance_ratio": lambda: ag.plot_variance_ratio(adata, n_pcs=30),
        "ggann.plot_clustermap": lambda: ag.plot_clustermap(adata, MARKERS, group_by=GROUP),
        "ggann.plot_upset": lambda: ag.plot_upset(marker_sets, min_cardinality=1),
    }


def _save(obj, path):
    """ggplot -> .save; marsilea/PyComplexHeatmap escape hatches -> current figure."""
    if isinstance(obj, p9.ggplot):
        obj.save(path, width=5.5, height=4, dpi=80, verbose=False)
    elif hasattr(obj, "save"):  # marsilea Upset
        obj.save(path, dpi=80)
    else:  # PyComplexHeatmap ClusterMapPlotter rendered onto the current figure
        plt.savefig(path, dpi=80, bbox_inches="tight")
    plt.close("all")


def main():
    adata = _adata()
    out = os.path.join(os.path.dirname(__file__), "..", "docs", "images", "api")
    os.makedirs(out, exist_ok=True)
    for name, builder in _examples(adata).items():
        _save(builder(), os.path.join(out, f"{name}.png"))
        print("wrote", name)


if __name__ == "__main__":
    main()
