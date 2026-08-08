"""Plots from differential-expression / marker results (``rank_genes_groups``).

scanpy stores marker results in ``adata.uns['rank_genes_groups']`` and exposes a
tidy view via ``sc.get.rank_genes_groups_df``. These helpers reuse that (never
parsing the recarrays by hand), pick the top markers, and delegate to the
existing ``plot_dotplot`` / ``plot_matrixplot`` -- or, for the volcano, to
plotnine-extra's ``ggvolcano``.
"""

from __future__ import annotations

import plotnine_extra as pe
from plotnine import (
    aes,
    geom_hline,
    geom_point,
    ggplot,
    labs,
    scale_color_manual,
    scale_x_log10,
)

from ._compat import renamed_keyword
from .plots import plot_dotplot, plot_matrixplot
from .theme import theme_ggann

# Conventional volcano colours: down = blue, non-significant = grey, up = red.
# Keys must match the categories ``plotnine_extra.ggvolcano`` maps to its colour
# aesthetic; unused keys are ignored, so a one-sided result (only "up") is fine.
_VOLCANO_COLORS = {"down": "#3B4CC0", "not significant": "#B8B8B8", "up": "#B40426"}

__all__ = [
    "rank_genes_df",
    "plot_rank_genes_dotplot",
    "plot_rank_genes_matrixplot",
    "plot_volcano",
    "plot_ma",
]

# MA-plot significance colours: significant = red, the rest grey.
_MA_COLORS = {True: "#B40426", False: "#B8B8B8"}


def _require_de(adata, key: str):
    if key not in adata.uns:
        raise KeyError(
            f"adata.uns[{key!r}] not found; run sc.tl.rank_genes_groups(adata, ...) first."
        )


def _de_group_by(adata, key: str, group_by):
    if group_by is not None:
        return group_by
    return adata.uns[key]["params"]["groupby"]


def rank_genes_df(adata, group=None, key: str = "rank_genes_groups", **kwargs):
    """Tidy DE table for ``group`` (or all groups), via ``sc.get.rank_genes_groups_df``.

    Columns: ``group, names, scores, logfoldchanges, pvals, pvals_adj``. Extra
    kwargs (``pval_cutoff``, ``log2fc_min``, ...) pass through to scanpy.
    """
    import scanpy as sc

    _require_de(adata, key)
    return sc.get.rank_genes_groups_df(adata, group=group, key=key, **kwargs)


def _top_genes(adata, n_genes: int, key: str) -> list[str]:
    df = rank_genes_df(adata, group=None, key=key)
    top = df.groupby("group", observed=True).head(n_genes)
    return list(dict.fromkeys(top["names"]))  # unique, group-ordered


def plot_rank_genes_dotplot(
    adata,
    n_genes: int = 5,
    group_by: str | None = None,
    key: str = "rank_genes_groups",
    **kwargs,
):
    """Dotplot of the top ``n_genes`` markers per group (``sc.pl.rank_genes_groups_dotplot``)."""
    _require_de(adata, key)
    group_by = renamed_keyword(
        group_by,
        kwargs.pop("groupby", None),
        name="group_by",
        legacy_name="groupby",
        default=None,
    )
    group_by = _de_group_by(adata, key, group_by)
    genes = _top_genes(adata, n_genes, key)
    return plot_dotplot(adata, genes, group_by, **kwargs)


def plot_rank_genes_matrixplot(
    adata,
    n_genes: int = 5,
    group_by: str | None = None,
    key: str = "rank_genes_groups",
    standard_scale: str | None = "var",
    **kwargs,
):
    """Matrixplot (group-mean tiles) of the top ``n_genes`` markers per group.

    Matches ``sc.pl.rank_genes_groups_matrixplot`` (group summary), not the
    per-cell ``..._heatmap``.
    """
    _require_de(adata, key)
    group_by = renamed_keyword(
        group_by,
        kwargs.pop("groupby", None),
        name="group_by",
        legacy_name="groupby",
        default=None,
    )
    group_by = _de_group_by(adata, key, group_by)
    genes = _top_genes(adata, n_genes, key)
    return plot_matrixplot(
        adata, genes, group_by, standard_scale=standard_scale, **kwargs
    )


def plot_volcano(
    adata,
    group: str,
    key: str = "rank_genes_groups",
    lfc: float = 1.0,
    padj: float = 0.05,
    label_top: int = 10,
    **kwargs,
):
    """Volcano plot (log2FC vs adjusted p-value) for one group's markers.

    Reuses plotnine-extra's ``ggvolcano``; ``lfc``/``padj`` set the fold-change /
    significance cutoffs and ``label_top`` labels the strongest genes. Requires a
    ``rank_genes_groups`` computed with a method that reports p-values and
    log-fold-changes (``wilcoxon`` / ``t-test``, not ``logreg``).
    """
    _require_de(adata, key)
    df = rank_genes_df(adata, group=group, key=key)
    missing = {"logfoldchanges", "pvals_adj"} - set(df.columns)
    if missing:
        raise ValueError(
            f"rank_genes_groups result is missing {sorted(missing)} needed for a volcano; "
            "re-run sc.tl.rank_genes_groups(adata, ..., method='wilcoxon' or 't-test')."
        )
    return (
        pe.ggvolcano(
            df,
            x="logfoldchanges",
            y="pvals_adj",
            label="names",
            p_cutoff=padj,
            fc_cutoff=lfc,
            label_top=label_top,
            **kwargs,
        )
        + scale_color_manual(values=_VOLCANO_COLORS)
        + labs(x="log2 fold change", y="-log10(adjusted p-value)")
        + theme_ggann()
    )


def plot_ma(
    data,
    *,
    mean: str = "baseMean",
    lfc: str = "log2FoldChange",
    pval: str = "padj",
    padj: float = 0.05,
    label: str | None = None,
    label_top: int = 0,
    size: float = 1.2,
):
    """MA plot of a (pseudobulk) differential-expression table.

    Plots mean expression on the x axis (log-scaled, the "A" of MA) against the
    log2 fold change on the y axis (the "M"), the standard diagnostic for a
    pseudobulk DE run (PyDESeq2 / decoupler / edgeR-style results). Genes with
    ``pval < padj`` are highlighted; the default column names follow PyDESeq2
    (``baseMean`` / ``log2FoldChange`` / ``padj``) but any table works via the
    ``mean`` / ``lfc`` / ``pval`` arguments.

    ``data`` is a :class:`pandas.DataFrame` (gene name in the index or a column
    named by ``label``). Set ``label_top=N`` to annotate the ``N`` genes with the
    largest absolute fold change among the significant ones.
    """
    missing = [c for c in (mean, lfc, pval) if c not in data.columns]
    if missing:
        raise KeyError(f"plot_ma: columns {missing} not in the results table.")

    # drop rows with no mean expression (can't sit on a log x axis); assign the
    # significance flag without a second full-frame copy
    df = data[data[mean] > 0].assign(
        significant=lambda d: (d[pval] < padj) & d[pval].notna()
    )

    plot = (
        ggplot(df, aes(mean, lfc, color="significant"))
        + geom_point(size=size, alpha=0.6, stroke=0)
        + geom_hline(yintercept=0, linetype="dashed", color="#4d4d4d")
        + scale_x_log10()
        + scale_color_manual(values=_MA_COLORS, labels={True: "sig.", False: "n.s."})
        + labs(x="mean expression", y="log2 fold change", color=f"{pval} < {padj}")
        + theme_ggann()
    )
    if label_top:
        sig = df[df["significant"]]
        top = sig.sort_values(lfc, key=lambda s: s.abs(), ascending=False).head(
            label_top
        )
        if label is not None:
            top = top.assign(_lab=top[label])
        else:
            top = top.assign(_lab=top.index.astype(str))
        plot = plot + pe.geom_text_repel(
            aes(mean, lfc, label="_lab"),
            data=top,
            size=8,
            color="black",
            inherit_aes=False,
        )
    return plot
