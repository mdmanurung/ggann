"""Standalone hierarchical dendrogram of groups (``sc.pl.dendrogram``).

scanpy stores a dendrogram under ``adata.uns['dendrogram_<groupby>']`` (computed
by ``sc.tl.dendrogram``), including the scipy plotting coordinates
(``icoord`` / ``dcoord``) and the leaf order (``ivl``). This module renders those
coordinates as an ordinary :class:`plotnine.ggplot` tree, so the dendrogram
composes with ggann's theme and can stand on its own -- unlike scanpy, where the
tree only appears bolted onto a dotplot / matrixplot.
"""

from __future__ import annotations

import pandas as pd
import plotnine_extra as pe
from plotnine import (
    aes,
    geom_line,
    ggplot,
    labs,
    scale_x_continuous,
    scale_y_continuous,
)

from .theme import theme_ggann

__all__ = ["plot_dendrogram"]


def _dendrogram_info(adata, group_by: str, key: str | None):
    key = key or f"dendrogram_{group_by}"
    if key not in adata.uns:
        import scanpy as sc

        result = sc.tl.dendrogram(adata, groupby=group_by, key_added=key, inplace=False)
        return result["dendrogram_info"]
    return adata.uns[key]["dendrogram_info"]


def plot_dendrogram(adata, group_by: str, *, key: str | None = None, orientation: str = "top"):
    """Hierarchical tree relating the categories of ``group_by`` (``sc.pl.dendrogram``).

    Reuses coordinates stored by ``sc.tl.dendrogram``. If they are absent, the
    coordinates are computed without modifying ``adata``. ``orientation='top'``
    draws leaves along the x axis; ``orientation='left'`` draws them along y.

    Parameters
    ----------
    adata
        Annotated data matrix.
    group_by : str
        Categorical observation column.
    key : str, optional
        ``adata.uns`` dendrogram key.
    orientation : {"top", "left"}
        Direction of the rendered tree.

    Returns
    -------
    plotnine.ggplot
        Composable dendrogram.

    Raises
    ------
    KeyError
        If grouping data or the requested stored result is missing.
    ValueError
        If ``orientation`` is unsupported or scanpy rejects grouping data.

    Notes
    -----
    Missing coordinates are computed on a temporary AnnData copy; input ownership
    is preserved.

    Examples
    --------
    >>> p = plot_dendrogram(adata, group_by="cell_type")
    """
    if orientation not in {"top", "left"}:
        raise ValueError("orientation must be 'top' or 'left'.")
    info = _dendrogram_info(adata, group_by, key)
    icoord, dcoord, ivl = info["icoord"], info["dcoord"], info["ivl"]

    # each link is a "|‾‾|" of 4 points; keep them as one grouped path
    rows = []
    for link, (xs, ys) in enumerate(zip(icoord, dcoord)):
        for x, y in zip(xs, ys):
            rows.append({"link": link, "pos": x, "height": y})
    seg = pd.DataFrame(rows)

    # scipy places the i-th leaf at x = 10*i + 5
    leaf_pos = [10 * i + 5 for i in range(len(ivl))]

    if orientation == "top":
        return (
            ggplot(seg, aes("pos", "height", group="link"))
            + geom_line()
            + scale_x_continuous(breaks=leaf_pos, labels=list(ivl))
            + labs(x="", y="distance")
            + theme_ggann()
            # rotate the leaf labels so long category names stay legible
            + pe.rotate_x_text(90)
        )
    # 'left' -- swap axes so leaves run down the y axis
    return (
        ggplot(seg, aes("height", "pos", group="link"))
        + geom_line()
        + scale_y_continuous(breaks=leaf_pos, labels=list(ivl))
        + labs(x="distance", y="")
        + theme_ggann()
    )
