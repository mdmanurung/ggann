"""Plots from differential-expression / marker results (``rank_genes_groups``).

scanpy stores marker results in ``adata.uns['rank_genes_groups']`` and exposes a
tidy view via ``sc.get.rank_genes_groups_df``. These helpers reuse that (never
parsing the recarrays by hand), pick the top markers, and delegate to the
existing ``plot_dotplot`` / ``plot_matrixplot``. Volcano and MA plots are built
from ordinary plotnine layers, so their prepared tables remain inspectable.
"""

from __future__ import annotations

import warnings
from typing import Sequence, cast

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

from ._compat import SKIP_GGANN_FRAMES, accepts_groupby
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
        If ``adata.uns[key]`` is absent or ``group`` is not a stored group
        (scanpy raises this ``KeyError``).
    TypeError
        If a filter such as ``pval_cutoff`` is not numeric (scanpy raises this).

    Notes
    -----
    Reads the stored result without modifying ``adata``. With a single group,
    scanpy drops the ``group`` column. A ``logreg`` result has only ``names``
    and ``scores``. ``pbmc3k_processed`` ships a result without
    ``logfoldchanges``; compute a fresh one into a separate key
    (``sc.tl.rank_genes_groups(adata, "louvain", method="wilcoxon",
    key_added="rank_wilcoxon")``) and pass ``key="rank_wilcoxon"``.

    Examples
    --------
    >>> table = rank_genes_df(adata, group="T cells")
    """
    import scanpy as sc

    _require_de(adata, key)
    return sc.get.rank_genes_groups_df(adata, group=group, key=key, **kwargs)


def _count(name: str, value, minimum: int) -> int:
    """Return ``value`` as an int >= ``minimum`` or raise a clear ValueError."""
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)):
        raise ValueError(f"{name} must be an integer, got {value!r}.")
    if value < minimum:
        raise ValueError(f"{name} must be at least {minimum}, got {value!r}.")
    return int(value)


def _positive_n_genes(n_genes) -> int:
    if isinstance(n_genes, (int, np.integer)) and not isinstance(n_genes, bool) and n_genes < 0:
        raise ValueError(
            f"n_genes must be a positive integer, got {n_genes!r}; scanpy's bottom-ranked "
            "selection (negative n_genes) is not supported."
        )
    return _count("n_genes", n_genes, 1)


def _marker_table(adata, n_genes, key, groups, min_logfoldchange) -> pd.DataFrame:
    """Top ``n_genes`` per marker group as ``[group, names]``, in group order."""
    n_genes = _positive_n_genes(n_genes)
    stored = list(adata.uns[key]["names"].dtype.names)
    if groups is None:
        groups = stored
    elif isinstance(groups, str):
        groups = [groups]
    unknown = [g for g in groups if g not in stored]
    if unknown:
        raise KeyError(f"groups {unknown} not in adata.uns[{key!r}]; stored groups: {stored}.")
    kwargs = {}
    if min_logfoldchange is not None:
        if "logfoldchanges" not in adata.uns[key]:
            raise ValueError(
                f"min_logfoldchange needs logfoldchanges, which adata.uns[{key!r}] lacks "
                "(e.g. method='logreg')."
            )
        kwargs["log2fc_min"] = min_logfoldchange
    parts = []
    for g in dict.fromkeys(groups):
        df = cast(pd.DataFrame, rank_genes_df(adata, group=g, key=key, **kwargs))
        df = df[df["names"].notna()].head(n_genes)
        parts.append(pd.DataFrame({"group": g, "names": df["names"].to_numpy()}))
    markers = pd.concat(parts, ignore_index=True)
    if markers.empty:
        raise ValueError("No markers left after applying groups/min_logfoldchange.")
    return markers


_DE_VALUES = ("scores", "logfoldchanges", "pvals", "pvals_adj", "log10_pvals", "log10_pvals_adj")


def _de_statistic(adata, values_to_plot: str, key: str) -> pd.Series:
    """DE statistic indexed by ``(group, gene)`` for every stored group."""
    if values_to_plot not in _DE_VALUES:
        raise ValueError(
            f"values_to_plot must be one of {list(_DE_VALUES)}, got {values_to_plot!r}."
        )
    column = values_to_plot.removeprefix("log10_")
    if column not in adata.uns[key]:
        raise ValueError(
            f"adata.uns[{key!r}] has no {column!r} (method='logreg' reports scores only)."
        )
    parts = []
    for g in adata.uns[key]["names"].dtype.names:
        df = cast(pd.DataFrame, rank_genes_df(adata, group=g, key=key))
        parts.append(df.assign(group=g)[["group", "names", column]])
    stat = cast(pd.Series, pd.concat(parts).set_index(["group", "names"])[column])
    if values_to_plot.startswith("log10_"):
        stat = pd.Series(-np.log10(stat.to_numpy(dtype=float)), index=stat.index)
    return stat


def _rank_plot(plot_fn, aes_name, adata, n_genes, group_by, key, groups, min_lfc, values, kwargs):
    """Select markers at the DE boundary, then build/adjust the expression plot."""
    _require_de(adata, key)
    group_by = _de_group_by(adata, key, group_by)
    markers = _marker_table(adata, n_genes, key, groups, min_lfc)
    repeated = bool(markers["names"].duplicated(keep=False).any())
    if (repeated or values is not None) and kwargs.get("backend", "plotnine") != "plotnine":
        raise ValueError(
            "repeated markers and values_to_plot are only supported with backend='plotnine'."
        )
    if values is not None and kwargs.get("standard_scale") is not None:
        raise ValueError(
            "standard_scale rescales expression means and cannot be combined with values_to_plot."
        )
    stat = _de_statistic(adata, values, key) if values is not None else None
    plot = plot_fn(adata, list(dict.fromkeys(markers["names"])), group_by, **kwargs)
    data = cast(pd.DataFrame, plot.data)
    if stat is not None:
        index = pd.MultiIndex.from_arrays([data[group_by].astype(str), data["feature"].astype(str)])
        if not index.isin(stat.index).all():
            raise ValueError(
                f"adata.uns[{key!r}] lacks statistics for some plotted genes or groups; "
                "re-run sc.tl.rank_genes_groups(..., n_genes=adata.n_vars) to store all genes."
            )
        data = data.assign(mean_expression=stat.reindex(index).to_numpy())
        title = str(values).replace("_", " ").replace("pvals", "p-value")
        title = "log fold change" if values == "logfoldchanges" else title
        plot = plot + labs(**{aes_name: title})
    if repeated:
        # One column per (marker group, gene) occurrence. Values come from one
        # per-gene computation, so repeated columns are exact copies.
        count = markers["names"].map(markers["names"].value_counts())
        labels = np.where(
            count > 1, markers["names"] + " (" + markers["group"] + ")", markers["names"]
        )
        data = pd.concat(
            [
                data[data["feature"] == gene].assign(feature=label, marker_group=group)
                for gene, label, group in zip(markers["names"], labels, markers["group"])
            ],
            ignore_index=True,
        )
        data["feature"] = pd.Categorical(
            data["feature"], categories=list(dict.fromkeys(labels)), ordered=True
        )
    plot.data = data
    return plot


@accepts_groupby
def plot_rank_genes_dotplot(
    adata,
    n_genes: int = 10,
    group_by: str | None = None,
    key: str = "rank_genes_groups",
    *,
    groups: str | Sequence[str] | None = None,
    min_logfoldchange: float | None = None,
    values_to_plot: str | None = None,
    **kwargs,
):
    """Dotplot of the top markers per group.

    Matches ``sc.pl.rank_genes_groups_dotplot``: ``n_genes=10`` markers per
    group, in DE order. A gene that is a top marker of several groups gets one
    column per marker group, labelled ``"GENE (group)"``; ``plot.data`` then
    carries a ``marker_group`` column. Scanpy's bottom-ranked selection
    (negative ``n_genes``) is not supported.

    Parameters
    ----------
    adata
        AnnData containing a scanpy rank-genes result.
    n_genes : int
        Markers selected per group; must be positive.
    group_by : str, optional
        Observation grouping column; inferred from the stored result by default.
    key : str
        Key in ``adata.uns``.
    groups : str or sequence of str, optional
        Marker groups to select from; all stored groups by default. The y axis
        still shows every category of ``group_by``.
    min_logfoldchange : float, optional
        Keep only markers with ``logfoldchanges`` above this value, applied
        before the top ``n_genes`` are taken.
    values_to_plot : str, optional
        One of ``scores``, ``logfoldchanges``, ``pvals``, ``pvals_adj``,
        ``log10_pvals``, ``log10_pvals_adj``. Colours dots by this DE statistic
        instead of mean expression (dot size stays the fraction expressing).
        Needs a result stored with ``n_genes=adata.n_vars``. Cannot be combined
        with ``standard_scale``.
    **kwargs
        Passed to :func:`ggann.plot_dotplot`.

    Returns
    -------
    plotnine.ggplot
        Composable marker dotplot.

    Raises
    ------
    KeyError
        If the stored result, a group, genes, or grouping column is missing.
    ValueError
        If ``n_genes`` is not a positive integer, the result lacks a requested
        statistic, or forwarded plotting options are invalid.

    Notes
    -----
    Only selected marker genes are projected for aggregation; input is unchanged.

    Examples
    --------
    >>> p = plot_rank_genes_dotplot(adata, n_genes=3)
    """
    return _rank_plot(
        plot_dotplot,
        "color",
        adata,
        n_genes,
        group_by,
        key,
        groups,
        min_logfoldchange,
        values_to_plot,
        kwargs,
    )


@accepts_groupby
def plot_rank_genes_matrixplot(
    adata,
    n_genes: int = 10,
    group_by: str | None = None,
    key: str = "rank_genes_groups",
    standard_scale: str | None = None,
    *,
    groups: str | Sequence[str] | None = None,
    min_logfoldchange: float | None = None,
    values_to_plot: str | None = None,
    **kwargs,
):
    """Matrixplot (group-mean tiles) of the top ``n_genes`` markers per group.

    Matches ``sc.pl.rank_genes_groups_matrixplot`` (group summary, raw group
    means by default), not the per-cell ``..._heatmap``. Repeated markers get
    one labelled column per marker group, as in
    :func:`ggann.plot_rank_genes_dotplot`.

    Parameters
    ----------
    adata
        AnnData containing a scanpy rank-genes result.
    n_genes : int
        Markers selected per group; must be positive.
    group_by : str, optional
        Observation grouping column; inferred by default.
    key : str
        Key in ``adata.uns``.
    standard_scale : str, optional
        Scaling forwarded to :func:`ggann.plot_matrixplot`; ``None`` (default)
        keeps raw group means, as scanpy does.
    groups, min_logfoldchange, values_to_plot
        As in :func:`ggann.plot_rank_genes_dotplot`; tiles show the DE
        statistic when ``values_to_plot`` is given.
    **kwargs
        Additional matrixplot arguments.

    Returns
    -------
    plotnine.ggplot
        Composable marker matrixplot.

    Raises
    ------
    KeyError
        If the stored result, a group, genes, or grouping column is missing.
    ValueError
        If ``n_genes`` is not a positive integer, the result lacks a requested
        statistic, or forwarded plotting options are invalid.

    Notes
    -----
    Only selected marker genes are projected for aggregation; input is unchanged.

    Examples
    --------
    >>> p = plot_rank_genes_matrixplot(adata, n_genes=3)
    """
    return _rank_plot(
        plot_matrixplot,
        "fill",
        adata,
        n_genes,
        group_by,
        key,
        groups,
        min_logfoldchange,
        values_to_plot,
        {**kwargs, "standard_scale": standard_scale},
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

    ``lfc`` and ``padj`` set the fold-change and significance cutoffs.
    ``label_top`` labels the most significant genes that pass both cutoffs.
    Requires a ``rank_genes_groups`` computed with a method that reports p-values
    and log-fold-changes (``wilcoxon`` / ``t-test``, not ``logreg``).

    Parameters
    ----------
    adata
        AnnData containing a scanpy rank-genes result.
    group : str
        One group to display. A list or ``None`` is rejected.
    key : str
        Key in ``adata.uns``.
    lfc : float
        Absolute log2-fold-change threshold.
    padj : float
        Adjusted-p-value threshold.
    label_top : int
        Number of strongest genes to label.
    **kwargs
        Passed to the point layer. Defaults are ``size=1.5``, ``alpha=0.8``, and
        ``stroke=0``. ``color``/``fill`` are rejected because colour encodes the
        direction; change it with ``scale_color_manual``. Extreme fold changes
        (for example a group with near-zero mean) stretch the x axis; zoom with
        ``+ coord_cartesian(xlim=(-5, 5))``.

    Returns
    -------
    plotnine.ggplot
        Composable volcano plot.

    Raises
    ------
    KeyError
        If the stored result or group is absent.
    ValueError
        If ``group`` is not a single str, p-values/fold changes are unavailable,
        or thresholds are invalid.

    Notes
    -----
    Reads only a stored rank-genes table and does not mutate ``adata``.

    Examples
    --------
    >>> p = plot_volcano(adata, group="T cells")
    """
    _require_de(adata, key)
    if not isinstance(group, str):
        raise ValueError(
            f"group must be one group name (str), got {group!r}; a volcano shows one "
            "group's markers, never several pooled."
        )
    if {"color", "colour", "fill"} & kwargs.keys():
        raise ValueError(
            "plot_volcano colours points by direction; override with "
            "+ scale_color_manual(...) instead of a color= keyword."
        )
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
    significant = (pvalues < padj) & (np.abs(fold_change) >= lfc)
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
    point_kwargs = {"size": 1.5, "alpha": 0.8, "stroke": 0, **kwargs}
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
    lfc_col: str = "log2FoldChange",
    pval: str = "padj",
    lfc: float = 1.0,
    padj: float = 0.05,
    label: str | None = None,
    label_top: int = 0,
    size: float = 1.2,
):
    """MA plot of a (pseudobulk) differential-expression table.

    Plots mean expression on the x axis (log-scaled, the "A" of MA) against the
    log2 fold change on the y axis (the "M"), the standard diagnostic for a
    pseudobulk DE run (PyDESeq2 / decoupler / edgeR-style results). The default
    column names follow PyDESeq2 (``baseMean`` / ``log2FoldChange`` / ``padj``);
    any table works via ``mean`` / ``lfc_col`` / ``pval``.

    ``lfc`` and ``padj`` are thresholds, as in :func:`ggann.plot_volcano`: a gene
    is significant when ``pval < padj`` and ``|log2 fold change| >= lfc``, and is
    coloured up (red) or down (blue). Pass ``lfc=0`` to colour on ``padj`` alone.

    ``data`` is a :class:`pandas.DataFrame` (gene name in the index or a column
    named by ``label``). Set ``label_top=N`` to annotate the ``N`` genes with the
    largest absolute fold change among the significant ones.

    Parameters
    ----------
    data : pandas.DataFrame
        Differential-expression table.
    mean, lfc_col, pval : str
        Mean-expression, log-fold-change, and adjusted-p-value column names.
    lfc : float
        Absolute log2-fold-change threshold.
    padj : float
        Adjusted-p-value threshold (strict ``<``).
    label : str, optional
        Gene-label column; the index is used by default.
    label_top : int
        Number of significant genes to label; must be a non-negative integer.
    size : float
        Point size.

    Returns
    -------
    plotnine.ggplot
        Composable MA plot.

    Raises
    ------
    TypeError
        If ``data`` is not a DataFrame.
    KeyError
        If a required column is absent.
    ValueError
        If thresholds or ``label_top`` are invalid, or no gene has a positive
        mean.

    Warns
    -----
    UserWarning
        Genes with a missing or non-positive ``mean`` cannot sit on the log x
        axis; they are dropped and the count is reported.

    Notes
    -----
    A prepared copy is attached to the plot; ``data`` is not modified.

    Examples
    --------
    >>> p = plot_ma(results)
    """
    if not isinstance(data, pd.DataFrame):
        raise TypeError(f"plot_ma needs a pandas DataFrame, got {type(data).__name__}.")
    missing = [c for c in (mean, lfc_col, pval, *([label] if label else [])) if c not in data]
    if missing:
        raise KeyError(f"plot_ma: columns {missing} not in the results table.")
    if isinstance(lfc, bool) or not np.isfinite(lfc) or lfc < 0:
        raise ValueError("lfc must be a finite non-negative number.")
    if isinstance(padj, bool) or not np.isfinite(padj) or not 0 < padj <= 1:
        raise ValueError("padj must be greater than zero and at most one.")
    label_top = _count("label_top", label_top, 0)

    # genes without a positive mean cannot sit on a log x axis: drop and say how many
    positive = (data[mean] > 0).to_numpy()
    if not positive.all():
        warnings.warn(
            f"plot_ma: dropped {int((~positive).sum())} of {len(data)} genes with missing "
            f"or non-positive {mean!r} (cannot be drawn on a log axis).",
            UserWarning,
            skip_file_prefixes=SKIP_GGANN_FRAMES,
        )
    if not positive.any():
        raise ValueError(f"plot_ma: no gene has a positive {mean!r}.")
    frame = data.loc[positive]
    fold_change = frame[lfc_col].to_numpy(dtype=float)
    significant = ((frame[pval] < padj) & frame[pval].notna()).to_numpy() & (
        np.abs(fold_change) >= lfc
    )
    direction = np.select(
        [significant & (fold_change < 0), significant & (fold_change > 0)],
        ["down", "up"],
        default="not significant",
    )
    df = frame.assign(
        significant=significant,
        _significance=pd.Categorical(direction, categories=list(_VOLCANO_COLORS), ordered=True),
    )

    plot = (
        ggplot(df, aes(mean, lfc_col, color="_significance"))
        + geom_point(size=size, alpha=0.8, stroke=0)
        + geom_hline(yintercept=0, linetype="dashed", color="#4d4d4d")
        + scale_x_log10()
        + scale_color_manual(values=_VOLCANO_COLORS, breaks=list(_VOLCANO_COLORS), drop=False)
        + labs(
            x="mean expression",
            y="log2 fold change",
            color=f"{pval} < {padj:g}\n|log2FC| ≥ {lfc:g}",
        )
        + _family_theme("standard")
    )
    if lfc:
        plot = plot + geom_hline(yintercept=(-lfc, lfc), linetype="dashed", color="#4d4d4d")
    if label_top:
        sig = df[df["significant"]]
        top = sig.sort_values(lfc_col, key=lambda s: s.abs(), ascending=False).head(label_top)
        if label is not None:
            top = top.assign(_lab=top[label])
        else:
            top = top.assign(_lab=top.index.astype(str))
        plot = plot + pe.geom_text_repel(
            aes(mean, lfc_col, label="_lab"),
            data=top,
            size=8,
            color="black",
            inherit_aes=False,
        )
    return plot
