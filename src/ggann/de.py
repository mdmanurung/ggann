"""Plots from differential-expression / marker results (``rank_genes_groups``).

scanpy stores marker results in ``adata.uns['rank_genes_groups']`` and exposes a
tidy view via ``sc.get.rank_genes_groups_df``. These helpers reuse that (never
parsing the recarrays by hand), pick the top markers, and delegate to the
existing ``plot_dotplot`` / ``plot_matrixplot``. Volcano and MA plots are built
from ordinary plotnine layers, so their prepared tables remain inspectable.
"""

from __future__ import annotations

from typing import cast

import numpy as np
import pandas as pd
import plotnine_extra as pe
from plotnine import (
    aes,
    geom_hline,
    geom_point,
    geom_vline,
    ggplot,
    labs,
    scale_color_manual,
    scale_x_log10,
)

from ._compat import renamed_keyword
from .plots import plot_dotplot, plot_matrixplot
from .publication import _family_theme

# Conventional volcano colours: down = blue, non-significant = grey, up = red.
# The ordered vocabulary remains stable even when a result is one-sided.
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

    Parameters
    ----------
    adata
        AnnData containing a scanpy rank-genes result.
    group : str or sequence of str, optional
        Group or groups to return; ``None`` returns all.
    key : str
        Key in ``adata.uns``.
    **kwargs
        Passed to ``scanpy.get.rank_genes_groups_df``.

    Returns
    -------
    pandas.DataFrame
        Tidy ranked-gene table.

    Raises
    ------
    KeyError
        If ``adata.uns[key]`` is absent.
    ValueError
        If scanpy rejects a group or filter.

    Notes
    -----
    Reads the stored result without modifying ``adata``.

    Examples
    --------
    >>> table = rank_genes_df(adata, group="T cells")
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
    """Dotplot of the top markers per group.

    Parameters
    ----------
    adata
        AnnData containing a scanpy rank-genes result.
    n_genes : int
        Markers selected per group.
    group_by : str, optional
        Observation grouping column; inferred from the stored result by default.
    key : str
        Key in ``adata.uns``.
    **kwargs
        Passed to :func:`ggann.plot_dotplot`.

    Returns
    -------
    plotnine.ggplot
        Composable marker dotplot.

    Raises
    ------
    KeyError
        If the stored result, genes, or grouping column is missing.
    ValueError
        If the result or forwarded plotting options are invalid.

    Notes
    -----
    Only selected marker genes are projected for aggregation; input is unchanged.

    Examples
    --------
    >>> p = plot_rank_genes_dotplot(adata, n_genes=3)
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

    Parameters
    ----------
    adata
        AnnData containing a scanpy rank-genes result.
    n_genes : int
        Markers selected per group.
    group_by : str, optional
        Observation grouping column; inferred by default.
    key : str
        Key in ``adata.uns``.
    standard_scale : str, optional
        Scaling forwarded to :func:`ggann.plot_matrixplot`.
    **kwargs
        Additional matrixplot arguments.

    Returns
    -------
    plotnine.ggplot
        Composable marker matrixplot.

    Raises
    ------
    KeyError
        If the stored result, genes, or grouping column is missing.
    ValueError
        If the result or forwarded plotting options are invalid.

    Notes
    -----
    Only selected marker genes are projected for aggregation; input is unchanged.

    Examples
    --------
    >>> p = plot_rank_genes_matrixplot(adata, n_genes=3)
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
    return plot_matrixplot(adata, genes, group_by, standard_scale=standard_scale, **kwargs)


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

    ``lfc`` and ``padj`` set the fold-change and significance cutoffs.
    ``label_top`` labels the most significant genes that pass both cutoffs.
    Requires a ``rank_genes_groups`` computed with a method that reports p-values
    and log-fold-changes (``wilcoxon`` / ``t-test``, not ``logreg``).

    Parameters
    ----------
    adata
        AnnData containing a scanpy rank-genes result.
    group : str
        Group to display.
    key : str
        Key in ``adata.uns``.
    lfc : float
        Absolute log2-fold-change threshold.
    padj : float
        Adjusted-p-value threshold.
    label_top : int
        Number of strongest genes to label.
    **kwargs
        Passed to the point layer. Defaults are ``size=1.2``, ``alpha=0.6``, and
        ``stroke=0``.

    Returns
    -------
    plotnine.ggplot
        Composable volcano plot.

    Raises
    ------
    KeyError
        If the stored result or group is absent.
    ValueError
        If p-values/fold changes are unavailable or thresholds are invalid.

    Notes
    -----
    Reads only a stored rank-genes table and does not mutate ``adata``.

    Examples
    --------
    >>> p = plot_volcano(adata, group="T cells")
    """
    _require_de(adata, key)
    df = cast(pd.DataFrame, rank_genes_df(adata, group=group, key=key))
    missing = {"logfoldchanges", "pvals_adj"} - set(df.columns)
    if missing:
        raise ValueError(
            f"rank_genes_groups result is missing {sorted(missing)} needed for a volcano; "
            "re-run sc.tl.rank_genes_groups(adata, ..., method='wilcoxon' or 't-test')."
        )
    if isinstance(lfc, bool) or not np.isfinite(lfc) or lfc < 0:
        raise ValueError("lfc must be a finite non-negative number.")
    if isinstance(padj, bool) or not np.isfinite(padj) or not 0 < padj <= 1:
        raise ValueError("padj must be greater than zero and at most one.")
    if isinstance(label_top, bool) or not isinstance(label_top, int) or label_top < 0:
        raise ValueError("label_top must be a non-negative integer.")
    if label_top and "names" not in df:
        raise ValueError("rank_genes_groups result is missing ['names'] needed for labels.")

    pvalues = np.asarray(cast(pd.Series, df.loc[:, "pvals_adj"]), dtype=float)
    positive = pvalues[np.isfinite(pvalues) & (pvalues > 0)]
    floor = max(
        np.finfo(float).tiny,
        float(np.min(positive)) / 10 if len(positive) else np.finfo(float).tiny,
    )
    fold_change = np.asarray(cast(pd.Series, df.loc[:, "logfoldchanges"]), dtype=float)
    significant = (pvalues <= padj) & (np.abs(fold_change) >= lfc)
    direction = np.select(
        [significant & (fold_change < 0), significant & (fold_change > 0)],
        ["down", "up"],
        default="not significant",
    )
    prepared = df.assign(
        _neg_log10_padj=-np.log10(np.clip(pvalues, floor, 1)),
        _significance=pd.Categorical(
            direction,
            categories=list(_VOLCANO_COLORS),
            ordered=True,
        ),
    )
    point_kwargs = {"size": 1.2, "alpha": 0.6, "stroke": 0, **kwargs}
    plot = (
        ggplot(prepared, aes("logfoldchanges", "_neg_log10_padj", color="_significance"))
        + geom_point(**point_kwargs)
        + geom_hline(yintercept=-np.log10(padj), linetype="dashed", color="#4D4D4D")
        + geom_vline(xintercept=(-lfc, lfc), linetype="dashed", color="#4D4D4D")
        + scale_color_manual(
            values=_VOLCANO_COLORS,
            breaks=list(_VOLCANO_COLORS),
            labels={
                "down": "down",
                "not significant": "not significant",
                "up": "up",
            },
            drop=False,
        )
        + labs(
            x="log2 fold change",
            y="-log10(adjusted p-value)",
            color=f"adjusted p ≤ {padj:g}\n|log2FC| ≥ {lfc:g}",
        )
        + _family_theme("standard")
    )
    if label_top:
        top = (
            prepared.loc[significant]
            .assign(_abs_lfc=lambda frame: frame["logfoldchanges"].abs())
            .sort_values(["pvals_adj", "_abs_lfc"], ascending=[True, False])
            .head(label_top)
        )
        plot = plot + pe.geom_text_repel(
            aes("logfoldchanges", "_neg_log10_padj", label="names"),
            data=top,
            inherit_aes=False,
            size=7,
            color="black",
        )
    return plot


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

    Parameters
    ----------
    data : pandas.DataFrame
        Differential-expression table.
    mean, lfc, pval : str
        Mean-expression, log-fold-change, and adjusted-p-value columns.
    padj : float
        Significance threshold.
    label : str, optional
        Gene-label column; the index is used by default.
    label_top : int
        Number of significant genes to label.
    size : float
        Point size.

    Returns
    -------
    plotnine.ggplot
        Composable MA plot.

    Raises
    ------
    KeyError
        If a required column is absent.
    ValueError
        If thresholds or label selection are invalid.

    Notes
    -----
    A prepared copy is attached to the plot; ``data`` is not modified.

    Examples
    --------
    >>> p = plot_ma(results)
    """
    missing = [c for c in (mean, lfc, pval) if c not in data.columns]
    if missing:
        raise KeyError(f"plot_ma: columns {missing} not in the results table.")

    # drop rows with no mean expression (can't sit on a log x axis); assign the
    # significance flag without a second full-frame copy
    df = data[data[mean] > 0].assign(significant=lambda d: (d[pval] < padj) & d[pval].notna())

    plot = (
        ggplot(df, aes(mean, lfc, color="significant"))
        + geom_point(size=size, alpha=0.6, stroke=0)
        + geom_hline(yintercept=0, linetype="dashed", color="#4d4d4d")
        + scale_x_log10()
        + scale_color_manual(values=_MA_COLORS, labels={True: "sig.", False: "n.s."})
        + labs(x="mean expression", y="log2 fold change", color=f"{pval} < {padj}")
        + _family_theme("standard")
    )
    if label_top:
        sig = df[df["significant"]]
        top = sig.sort_values(lfc, key=lambda s: s.abs(), ascending=False).head(label_top)
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
