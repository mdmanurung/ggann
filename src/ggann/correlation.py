"""Group/sample correlation heatmap from pseudobulk mean-expression profiles.

Inspired by DOTools' ``correlation`` and marsilea's annotated heatmaps: aggregate
mean expression per group (a pseudobulk profile per cell type / sample / batch),
correlate those profiles, and draw the group-by-group correlation as a plotnine
tile plot. Optional hierarchical clustering reorders the axes so related groups
sit together (marsilea/ComplexHeatmap style), and values can be annotated in-cell.
"""

from __future__ import annotations

import warnings
from typing import Sequence

import numpy as np
import pandas as pd
import plotnine_extra as pe
from plotnine import (
    aes,
    coord_equal,
    element_blank,
    geom_tile,
    ggplot,
    labs,
    scale_fill_cmap,
    scale_fill_gradient2,
    theme,
)

from ._aggregate import group_means
from ._annotation import annotation_threshold, geom_contrast_text
from ._expression import ordered_unique
from .publication import _active_style, _family_theme

__all__ = ["plot_correlation"]


def _default_genes(adata) -> list[str]:
    """Highly-variable genes when flagged, else every var (selection, not extraction)."""
    if "highly_variable" in adata.var.columns:
        hv = adata.var_names[np.asarray(adata.var["highly_variable"], dtype=bool)]
        if len(hv) > 1:
            return list(hv)
    genes = list(adata.var_names)
    warnings.warn(
        f"plot_correlation: no 'highly_variable' genes flagged; correlating all "
        f"{len(genes)} genes. Pass genes=... to restrict (faster on large data).",
        stacklevel=3,
    )
    return genes


def _cluster_order(corr: pd.DataFrame) -> list[str]:
    """Leaf order from average-linkage clustering on ``1 - correlation`` distance.

    A group with a constant (zero-variance) profile yields NaN correlations; those
    become the maximum distance (2.0) so the group is simply placed as an outlier
    instead of crashing scipy's ``linkage`` on non-finite values.
    """
    from scipy.cluster.hierarchy import leaves_list, linkage
    from scipy.spatial.distance import squareform

    dist = 1.0 - corr.to_numpy()
    dist = np.nan_to_num(dist, nan=2.0)  # NaN corr -> max distance, never non-finite
    np.fill_diagonal(dist, 0.0)
    dist = (dist + dist.T) / 2.0  # enforce symmetry for squareform
    order = leaves_list(linkage(squareform(dist, checks=False), method="average"))
    return [str(corr.columns[int(i)]) for i in order]


def _fill_scale(values: pd.Series, cmap: str | None):
    """Honest fill scale for correlations: diverging-at-zero when signs mix, else sequential.

    A diverging map stretched over an all-positive range would put its neutral
    colour mid-range (making a strong 0.7 look anti-correlated); so all-positive
    data gets a sequential map over its own range instead.
    """
    finite = values.to_numpy()
    finite = finite[np.isfinite(finite)]
    lo = float(finite.min()) if finite.size else -1.0
    if cmap is not None:
        style = _active_style()
        kwargs = {"na_value": style.missing_color} if style is not None else {}
        return scale_fill_cmap(cmap_name=cmap, limits=(lo, 1.0), **kwargs)
    style = _active_style()
    if lo < 0.0:  # mixed sign -> diverging, white centred on zero
        m = float(np.abs(finite).max()) if finite.size else 1.0
        if style is not None:
            return scale_fill_gradient2(
                low=style.diverging[0],
                mid=style.diverging[1],
                high=style.diverging[2],
                midpoint=0,
                limits=(-m, m),
                na_value=style.missing_color,
            )
        return scale_fill_cmap(cmap_name="RdBu_r", limits=(-m, m))
    if style is not None:
        return scale_fill_cmap(
            cmap_name=style.sequential_cmap,
            limits=(lo, 1.0),
            na_value=style.missing_color,
        )
    return scale_fill_cmap(cmap_name="YlOrRd", limits=(lo, 1.0))  # all-positive -> sequential


def plot_correlation(
    adata,
    group_by: str,
    *,
    genes: Sequence[str] | None = None,
    layer: str | None = None,
    use_raw: bool | None = None,
    method: str = "pearson",
    cluster: bool = True,
    annotate: bool | str = False,
    cmap: str | None = None,
    rasterized: bool = False,
):
    """Correlation heatmap between the mean-expression profiles of each group.

    ``genes`` defaults to the highly-variable set (or all genes, with a warning).
    ``method`` is any :meth:`pandas.DataFrame.corr` method (``"pearson"`` /
    ``"spearman"`` / ``"kendall"``). With ``cluster=True`` the axes are reordered by
    hierarchical clustering; ``annotate=True`` prints each correlation in its cell.
    ``cmap=None`` uses a sequential scale when all correlations are positive and a
    zero-centred diverging scale when signs differ. Pass a colormap name to set it
    explicitly.

    Parameters
    ----------
    adata
        Annotated data matrix.
    group_by : str
        Observation column defining mean profiles.
    genes : sequence of str, optional
        Genes included in profile correlations.
    layer, use_raw : optional
        Mutually exclusive expression source.
    method : {"pearson", "spearman", "kendall"}
        Pandas correlation method.
    cluster : bool
        Hierarchically order groups.
    annotate : bool or {"auto", "force"}
        Print correlation coefficients. ``"auto"`` labels only cells at least
        12 points wide and high; ``"force"`` always labels. ``True`` is an
        alias for ``"force"``.
    cmap : str, optional
        Explicit Matplotlib colormap.
    rasterized : bool, default=False
        Rasterize the heatmap tiles while retaining vector labels and guides.

    Returns
    -------
    plotnine.ggplot
        Composable correlation heatmap.

    Raises
    ------
    KeyError
        If genes, grouping metadata, or a layer are missing.
    ValueError
        If fewer than two genes resolve or method/source selection is invalid.

    Notes
    -----
    Requested genes are projected and sparse-aggregated before the small group
    correlation matrix is materialized. ``adata`` is unchanged.

    Examples
    --------
    >>> p = plot_correlation(adata, group_by="cell_type", genes=["CD3D", "NKG7"])
    """
    if genes is None:
        genes = _default_genes(adata)
    genes = ordered_unique(genes)
    if len(genes) < 2:
        raise ValueError("plot_correlation needs at least two genes to correlate profiles.")

    profiles = group_means(adata, genes, group_by, layer=layer, use_raw=use_raw).T  # genes x groups
    corr = profiles.corr(method=method)

    order = _cluster_order(corr) if cluster and corr.shape[0] >= 2 else list(corr.columns)
    corr = corr.loc[order, order]

    long = (
        corr.rename_axis(index="row", columns="col")
        .reset_index()
        .melt(id_vars="row", var_name="col", value_name="corr")
    )
    long["row"] = pd.Categorical(long["row"], categories=order, ordered=True)
    long["col"] = pd.Categorical(long["col"], categories=order, ordered=True)

    plot = (
        ggplot(long, aes("col", "row", fill="corr"))
        + geom_tile(raster=rasterized)
        + _fill_scale(long["corr"], cmap)
        + coord_equal()
        + labs(x="", y="", fill=f"{method}\ncorrelation")
        + _family_theme("matrix")
        + theme(axis_ticks=element_blank())
        + pe.rotate_x_text(45)
    )
    threshold = annotation_threshold(annotate)
    if threshold is not None:
        style = _active_style()
        plot = plot + geom_contrast_text(
            aes(label="corr"),
            format_string="{:.2f}",
            min_cell_pt=threshold,
            size=style.axis_text_size if style is not None else 7,
        )
    return plot
