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
    geom_path,
    ggplot,
    labs,
    scale_x_continuous,
    scale_y_continuous,
)

from ._compat import accepts_groupby
from ._grouping import _drop_missing_groups
from .publication import _family_theme

__all__ = ["plot_dendrogram"]


def _dendrogram_info(adata, group_by: str, key: str | None):
    key = key or f"dendrogram_{group_by}"
    if key not in adata.uns:
        import scanpy as sc

        if group_by not in adata.obs:
            raise KeyError(f"{group_by!r} is not a column of adata.obs.")
        column = adata.obs[group_by]
        if not isinstance(column.dtype, pd.CategoricalDtype):
            raise ValueError(
                f"{group_by!r} must be a categorical obs column; "
                f"convert it with adata.obs[{group_by!r}].astype('category')."
            )
        observed = column.cat.remove_unused_categories()
        if len(observed.cat.categories) < 2:
            raise ValueError(
                f"{group_by!r} has {len(observed.cat.categories)} observed categories; "
                "a dendrogram needs at least 2."
            )
        if len(observed.cat.categories) < len(column.cat.categories):
            # scanpy fails on unused categories; compute on a copy without them
            adata = adata.copy()
            adata.obs[group_by] = observed

        result = sc.tl.dendrogram(adata, groupby=group_by, key_added=key, inplace=False)
        if result is None:
            raise RuntimeError("scanpy did not return dendrogram coordinates.")
        return result["dendrogram_info"]
    return adata.uns[key]["dendrogram_info"]


@accepts_groupby
def plot_dendrogram(adata, group_by: str, *, key: str | None = None, orientation: str = "top"):
    """Hierarchical tree relating the categories of ``group_by`` (``sc.pl.dendrogram``).

    Reuses coordinates stored by ``sc.tl.dendrogram``. If they are absent, the
    coordinates are computed without modifying ``adata``; categories with no
    cells are ignored. A ``key`` that is not in ``adata.uns`` is not an error: the
    tree is recomputed and the key is not stored. ``orientation='top'`` draws
    leaves along the x axis; ``orientation='left'`` draws them along y.

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
        If ``group_by`` is not an ``obs`` column and no stored result is used.
    ValueError
        If ``orientation`` is unsupported, ``group_by`` is not categorical, or
        fewer than two categories have cells.

    Notes
    -----
    Missing coordinates are computed on a temporary AnnData copy; input ownership
    is preserved.

    Examples
    --------
    >>> p = plot_dendrogram(adata, group_by="cell_type")
    """
    adata = _drop_missing_groups(adata, group_by)
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
            + geom_path()
            + scale_x_continuous(breaks=leaf_pos, labels=list(ivl))
            + labs(x="", y="distance")
            + _family_theme("matrix")
            # rotate the leaf labels so long category names stay legible
            + pe.rotate_x_text(90)
        )
    # 'left' -- swap axes so leaves run down the y axis
    return (
        ggplot(seg, aes("height", "pos", group="link"))
        + geom_path()
        + scale_y_continuous(breaks=leaf_pos, labels=list(ivl))
        + labs(x="distance", y="")
        + _family_theme("matrix")
    )
