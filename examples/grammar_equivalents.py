"""Reference plotnine constructions for selected ggann helpers.

Each ``*_grammar`` function exposes the data preparation and plotnine layers
used for one helper-shaped figure. Helpers may also apply validation, ordering,
and styling that is omitted from these compact constructions.

``gganndata(adata, aes(...))`` does the name resolution and returns a real
``plotnine.ggplot`` whose ``.data`` is the tidy frame the layers draw from; for
grouped summaries the annplyr accessor (``adata.ap.summarize``) does the
aggregation. From there it is ordinary plotnine (``+ geom_* + scale_* + theme``).

``HELPERS`` maps each name to the corresponding convenience call so ``main()``
can render both and save them side by side under
``docs/images/compare/``. Escape hatches (``plot_clustermap`` via PyComplexHeatmap,
``plot_upset`` via marsilea) are not plotnine and have no grammar equivalent.

Run: ``python examples/grammar_equivalents.py``.
"""

from __future__ import annotations

import os

import matplotlib

matplotlib.use("Agg")

import annplyr as ap
import numpy as np
import pandas as pd
import plotnine_extra as pe
from plotnine import (
    aes,
    geom_boxplot,
    geom_col,
    geom_errorbar,
    geom_jitter,
    geom_line,
    geom_point,
    geom_tile,
    element_blank,
    geom_violin,
    ggplot,
    guide_legend,
    guides,
    labs,
    theme,
    scale_color_cmap,
    scale_fill_cmap,
    scale_size,
)

import ggann as ag
from ggann import gganndata

GROUP = "bulk_labels"
MARKERS = ["CD3D", "NKG7", "CST3", "GNLY"]


# --- plot_embedding(color=..., label=True) --------------------------------- #
def embedding_grammar(adata, color=GROUP, basis="umap", label=True):
    coords = ag.embedding_coords(adata, basis)
    xcol, ycol = coords.columns[:2]
    base = gganndata(adata, ag.aes(xcol, ycol, color=color))  # resolves the obs column
    plot = (
        base
        + geom_point(size=1.5, alpha=0.9)
        + ag.scale_color_obs(adata, color)
        + guides(color=guide_legend(override_aes={"size": 4}))
        + ag.theme_ggann()
        + theme(
            axis_text=element_blank(), axis_ticks=element_blank()
        )  # arbitrary UMAP units
    )
    if label:
        cents = (
            base.data.groupby(color, observed=True)[[xcol, ycol]].median().reset_index()
        )
        cents = cents.rename(columns={color: "label"})
        plot = plot + pe.geom_label_repel(
            aes(xcol, ycol, label="label"),
            data=cents,
            fill="white",
            color="black",
            inherit_aes=False,
        )
    return plot


# --- plot_features(features) ----------------------------------------------- #
def features_grammar(adata, features=("CD3D", "NKG7", "CST3", "GNLY"), basis="umap"):
    coords = ag.embedding_coords(adata, basis)
    xcol, ycol = coords.columns[:2]
    feats = list(features)
    wide = gganndata(adata, ag.aes(xcol, ycol, color=ag.gene(feats[0]))).data[
        [xcol, ycol]
    ]
    for f in feats:
        wide[f] = gganndata(adata, ag.aes(xcol, ycol, color=ag.gene(f))).data[f]
    long = wide.melt(
        id_vars=[xcol, ycol],
        value_vars=feats,
        var_name="feature",
        value_name="expression",
    )
    long["feature"] = pd.Categorical(long["feature"], categories=feats, ordered=True)
    long = long.sort_values("expression")
    return (
        ggplot(long, aes(xcol, ycol, color="expression"))
        + geom_point(size=1.2, alpha=0.9)
        + pe.facet_wrap("~feature")
        + scale_color_cmap(cmap_name="magma")
        + ag.theme_ggann()
    )


# --- plot_density(feature) ------------------------------------------------- #
def density_grammar(adata, feature="CD3D", basis="umap"):
    from pynebulosa import calculate_density

    coords = ag.embedding_coords(adata, basis)
    xcol, ycol = coords.columns[:2]
    resolved = gganndata(adata, ag.aes(xcol, ycol, color=ag.gene(feature))).data
    density = calculate_density(
        resolved[feature].to_numpy(float), resolved[[xcol, ycol]].to_numpy(float)
    )
    df = resolved.assign(density=density).sort_values("density")
    return (
        ggplot(df, aes(xcol, ycol, color="density"))
        + geom_point(size=1.5, alpha=0.9)
        + scale_color_cmap(cmap_name="magma")
        + labs(color="density")
        + ag.theme_ggann()
    )


def _order_group(df, group, adata):
    """Order a group column by its obs categorical order, as the helpers do."""
    col = adata.obs[group]
    cats = (
        list(col.cat.categories)
        if hasattr(col, "cat")
        else sorted(pd.unique(df[group]))
    )
    df = df.copy()
    df[group] = pd.Categorical(df[group], categories=cats, ordered=False)
    return df


def _grouped_means(adata, genes, group, cutoff=0.0):
    """Group-by-gene mean expression and fraction expressing, via annplyr."""
    mean = adata.ap.summarize(raw={g: ap.mean(ap.col(g)) for g in genes}, by=group)
    frac = adata.ap.summarize(
        raw={g: ap.mean(ap.col(g) > cutoff) for g in genes}, by=group
    )
    mean = mean.set_index(group)[genes].astype(float)
    frac = frac.set_index(group)[genes].astype(float)
    long = (
        mean.reset_index()
        .melt(id_vars=group, var_name="feature", value_name="mean")
        .merge(
            frac.reset_index().melt(
                id_vars=group, var_name="feature", value_name="frac"
            ),
            on=[group, "feature"],
        )
    )
    long["feature"] = pd.Categorical(
        long["feature"], categories=list(genes), ordered=True
    )
    return _order_group(long, group, adata)


# --- plot_dotplot(genes, group) -------------------------------------------- #
def dotplot_grammar(adata, genes=MARKERS, group=GROUP):
    long = _grouped_means(adata, list(genes), group)
    return (
        ggplot(long, aes("feature", group))
        + geom_point(aes(size="frac", color="mean"))
        + scale_color_cmap(cmap_name="Reds")
        + scale_size(range=(0.5, 8.0), labels=lambda xs: [f"{x:.0%}" for x in xs])
        + labs(x="", y="", color="mean\nexpression", size="fraction\nexpressing")
        + ag.theme_ggann()
        + pe.rotate_x_text(45)
    )


# --- plot_matrixplot(genes, group) ----------------------------------------- #
def matrixplot_grammar(adata, genes=MARKERS, group=GROUP):
    long = _grouped_means(adata, list(genes), group)
    return (
        ggplot(long, aes("feature", group, fill="mean"))
        + geom_tile()
        + scale_fill_cmap(cmap_name="viridis")
        + labs(x="", y="", fill="mean\nexpression")
        + ag.theme_ggann()
        + pe.rotate_x_text(45)
    )


# --- plot_violin(genes, group) --------------------------------------------- #
def violin_grammar(adata, gene="CD3D", group=GROUP):
    # one gene shown; for a facet-per-gene grid melt the genes and add facet_wrap
    d = _order_group(
        gganndata(adata, ag.aes(group, ag.gene(gene), fill=group)).data, group, adata
    )
    return (
        ggplot(d, aes(group, gene, fill=group))
        + geom_violin(scale="width")
        + geom_boxplot(width=0.12, fill="white", outlier_alpha=0.0, show_legend=False)
        + ag.scale_fill_obs(adata, group)
        + labs(x="", y="expression")
        + ag.theme_ggann()
        + pe.rotate_x_text(45)
    )


# --- plot_box(gene, group) ------------------------------------------------- #
def box_grammar(adata, gene="CD3D", group=GROUP):
    d = _order_group(
        gganndata(adata, ag.aes(group, ag.gene(gene), fill=group)).data, group, adata
    )
    return (
        ggplot(d, aes(group, gene, fill=group))
        + geom_boxplot(width=0.7, outlier_alpha=0.0)
        + geom_jitter(width=0.2, height=0.0, size=0.35, alpha=0.25, stroke=0)
        + ag.scale_fill_obs(adata, group)
        + labs(x="", y="expression")
        + ag.theme_ggann()
        + pe.rotate_x_text(45)
    )


# --- plot_expression_bar(gene, group) -------------------------------------- #
def bar_grammar(adata, gene="CD3D", group=GROUP):
    d = gganndata(adata, ag.aes(group, ag.gene(gene))).data
    s = (
        d.groupby(group, observed=True)[gene]
        .agg(mean="mean", sd="std", n="count")
        .reset_index()
    )
    s["se"] = s["sd"] / np.sqrt(s["n"])
    s["ymin"], s["ymax"] = s["mean"] - s["se"], s["mean"] + s["se"]
    s = _order_group(s, group, adata)
    return (
        ggplot(s, aes(group, "mean", fill=group))
        + geom_col(width=0.7)
        + geom_errorbar(aes(ymin="ymin", ymax="ymax"), width=0.3)
        + ag.scale_fill_obs(adata, group)
        + labs(x="", y="mean expression")
        + ag.theme_ggann()
        + pe.rotate_x_text(45)
    )


# --- plot_expression_line(gene, x, group) ---------------------------------- #
def line_grammar(adata, gene="CD3D", x="phase", group=GROUP):
    d = gganndata(adata, ag.aes(x, ag.gene(gene), color=group)).data
    s = d.groupby([x, group], observed=True)[gene].mean().reset_index(name="mean")
    return (
        ggplot(s, aes(x, "mean", color=group, group=group))
        + geom_line()
        + geom_point(size=2)
        + ag.scale_color_obs(adata, group)
        + labs(x=x, y="mean expression")
        + ag.theme_ggann()
    )


# --- plot_proportions(group, split_by) ------------------------------------- #
def proportions_grammar(adata, group=GROUP, split="phase"):
    d = gganndata(adata, ag.aes(split, fill=group)).data
    counts = d.groupby([split, group], observed=True).size().reset_index(name="n")
    counts["frac"] = counts.groupby(split, observed=True)["n"].transform(
        lambda x: x / x.sum()
    )
    counts = _order_group(counts, group, adata)
    return (
        ggplot(counts, aes(split, "frac", fill=group))
        + geom_col()
        + ag.scale_fill_obs(adata, group)
        + labs(x="", y="proportion")
        + ag.theme_ggann()
    )


# --- plot_correlation(group) ----------------------------------------------- #
def correlation_grammar(adata, group=GROUP, n_genes=None):
    genes = list(adata.var_names[np.asarray(adata.var["highly_variable"], bool)])
    if n_genes is not None:
        genes = genes[:n_genes]
    # Match plot_correlation's default source: raw when present, otherwise X.
    source = "raw" if adata.raw is not None else "x"
    means = adata.ap.summarize(
        **{source: {g: ap.mean(ap.col(g)) for g in genes}},
        by=group,
    )
    profiles = (
        means.set_index(group)[genes].apply(lambda c: np.asarray(c, float)).T
    )  # genes x groups
    corr = profiles.corr()
    long = (
        corr.rename_axis("row")
        .reset_index()
        .melt(id_vars="row", var_name="col", value_name="corr")
    )
    return (
        ggplot(long, aes("col", "row", fill="corr"))
        + geom_tile()
        + scale_fill_cmap(cmap_name="RdBu_r")
        + labs(x="", y="", fill="pearson\ncorrelation")
        + ag.theme_ggann()
        + pe.rotate_x_text(45)
    )


# Helper name -> reference construction used by the docs and smoke tests.
EQUIVALENTS = {
    "plot_embedding": embedding_grammar,
    "plot_features": features_grammar,
    "plot_density": density_grammar,
    "plot_dotplot": dotplot_grammar,
    "plot_matrixplot": matrixplot_grammar,
    "plot_violin": violin_grammar,
    "plot_box": box_grammar,
    "plot_expression_bar": bar_grammar,
    "plot_expression_line": line_grammar,
    "plot_proportions": proportions_grammar,
    "plot_correlation": correlation_grammar,
}

# Helper name -> corresponding convenience call.
HELPERS = {
    "plot_embedding": lambda ad: ag.plot_embedding(ad, "umap", color=GROUP, label=True),
    "plot_features": lambda ad: ag.plot_features(
        ad, ["CD3D", "NKG7", "CST3", "GNLY"], basis="umap"
    ),
    "plot_density": lambda ad: ag.plot_density(ad, "CD3D", basis="umap"),
    "plot_dotplot": lambda ad: ag.plot_dotplot(ad, MARKERS, GROUP),
    "plot_matrixplot": lambda ad: ag.plot_matrixplot(ad, MARKERS, GROUP),
    "plot_violin": lambda ad: ag.plot_violin(ad, ["CD3D"], GROUP),
    "plot_box": lambda ad: ag.plot_box(ad, ["CD3D"], GROUP),
    "plot_expression_bar": lambda ad: ag.plot_expression_bar(ad, ["CD3D"], GROUP),
    "plot_expression_line": lambda ad: ag.plot_expression_line(
        ad, ["CD3D"], x="phase", group_by=GROUP
    ),
    "plot_proportions": lambda ad: ag.plot_proportions(
        ad, GROUP, split_by="phase", position="fill"
    ),
    "plot_correlation": lambda ad: ag.plot_correlation(ad, GROUP),
}


def main():
    import scanpy as sc

    adata = sc.datasets.pbmc68k_reduced()
    out = os.path.join(os.path.dirname(__file__), "..", "docs", "images", "compare")
    os.makedirs(out, exist_ok=True)
    for name, grammar in EQUIVALENTS.items():
        for kind, plot in (
            ("helper", HELPERS[name](adata)),
            ("grammar", grammar(adata)),
        ):
            plot.save(
                os.path.join(out, f"{name}_{kind}.png"),
                width=5,
                height=4,
                dpi=80,
                verbose=False,
            )
        print("wrote", name)


if __name__ == "__main__":
    main()
