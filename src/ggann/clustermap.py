"""Clustered-heatmap escape hatch, delegating to PyComplexHeatmap.

Clustered heatmaps with dendrograms and stacked annotation bars are a grid-based
paradigm that does not fit the grammar of graphics, so they live *outside* the
``gganndata() + geom_*`` path on purpose. This one function is the bridge: it
builds an annplyr-tidied matrix and hands it to
:class:`PyComplexHeatmap.ClusterMapPlotter`, returning that plotter object.

``PyComplexHeatmap`` is an optional dependency; it is imported lazily so that
``import ggann`` never requires it.
"""

from __future__ import annotations

from typing import Iterable, Sequence

import pandas as pd

from ._aggregate import _standardize, group_means
from ._expression import densify_frame as _densify
from ._expression import ordered_unique, project_expression, resolve_source

__all__ = ["plot_clustermap"]


def _require_pch():
    try:
        import PyComplexHeatmap as pch
    except ImportError as exc:  # pragma: no cover - env-dependent
        raise ImportError(
            "plot_clustermap requires PyComplexHeatmap. "
            "Install it with `pip install ggann[heatmap]` or `pip install PyComplexHeatmap`."
        ) from exc
    return pch


def plot_clustermap(
    adata,
    genes: Sequence[str],
    group_by: str | None = None,
    *,
    layer: str | None = None,
    use_raw: bool | None = None,
    annotations: Iterable[str] | None = None,
    standard_scale: str | None = None,
    z_score: int | None = None,
    row_cluster: bool = True,
    col_cluster: bool = True,
    cmap: str = "viridis",
    show_rownames: bool = True,
    show_colnames: bool = True,
    plot: bool = True,
    **kwargs,
):
    """Clustered heatmap of ``genes`` (rows) across groups or cells (columns).

    * ``group_by`` given -> columns are aggregated group means (like
      ``sc.pl.heatmap`` on pseudobulk). Column annotations default to the group
      identity.
    * ``group_by=None`` -> columns are individual cells; ``annotations`` selects
      ``obs`` columns to show as top annotation bars.

    ``standard_scale`` (ggann-side, per gene/group) and ``z_score``
    (PyComplexHeatmap-side row/column z-score) are mutually exclusive -- applying
    both would normalize twice. Returns the
    :class:`PyComplexHeatmap.ClusterMapPlotter` instance.

    Parameters
    ----------
    adata
        Annotated data matrix.
    genes : sequence of str
        Genes displayed as rows.
    group_by : str, optional
        Observation column to aggregate into columns; otherwise show cells.
    layer, use_raw : optional
        Mutually exclusive expression source.
    annotations : iterable of str, optional
        Observation annotations for per-cell columns.
    standard_scale : {None, "var", "group", "zscore"}
        ggann-side scaling before clustering.
    z_score : {0, 1}, optional
        Backend row or column z-score.
    row_cluster, col_cluster : bool
        Enable hierarchical clustering.
    cmap : str
        Matplotlib colormap name.
    show_rownames, show_colnames : bool
        Show matrix labels.
    plot : bool
        Render immediately.
    **kwargs
        Passed to ``PyComplexHeatmap.ClusterMapPlotter``.

    Returns
    -------
    PyComplexHeatmap.ClusterMapPlotter
        Grid-backend plotter; this is not a plotnine object.

    Raises
    ------
    ImportError
        If the ``heatmap`` extra is unavailable.
    KeyError
        If requested genes, metadata, or a layer are missing.
    ValueError
        If sources are incompatible or both scaling mechanisms are selected.

    Notes
    -----
    Requested genes are projected before aggregation. Per-cell mode materializes a
    genes-by-cells table for clustering; ``adata`` remains unchanged.

    Examples
    --------
    >>> cm = plot_clustermap(adata, ["CD3D", "NKG7"], group_by="cell_type")
    """
    if standard_scale is not None and z_score is not None:
        raise ValueError(
            "standard_scale and z_score are mutually exclusive (both normalize the "
            "matrix); set at most one."
        )
    pch = _require_pch()
    genes = ordered_unique(genes)

    if group_by is not None:
        means = group_means(adata, genes, group_by, layer=layer, use_raw=use_raw)
        matrix = _standardize(means, standard_scale).T.loc[genes]
        ann_df = pd.DataFrame({group_by: matrix.columns}, index=matrix.columns)
    else:
        kind, lyr = resolve_source(adata, layer, use_raw)
        ann_cols = list(annotations) if annotations else []
        projected, genes = project_expression(adata, genes, kind=kind, layer=lyr, obs=ann_cols)
        wide = _densify(projected.ap.to_df(x=genes))[genes]
        wide.index = adata.obs_names.copy()
        wide = _standardize(wide, standard_scale)
        matrix = wide.T  # genes x cells
        ann_df = projected.ap.to_df(obs=ann_cols) if ann_cols else None

    top_annotation = None
    if ann_df is not None and not ann_df.empty:
        # pandas>=3 gives string columns the `str` dtype, which PyComplexHeatmap
        # does not recognise as categorical and then refuses to auto-pick a cmap.
        # Coerce non-numeric annotations to `category` so it colours them discretely.
        # `is_numeric_dtype` is True for bool, so coerce bool columns too -- they
        # are categorical (e.g. is_doublet), not a continuous scale.
        ann_df = ann_df.copy()
        for col in ann_df.columns:
            is_numeric = pd.api.types.is_numeric_dtype(ann_df[col])
            if not is_numeric or pd.api.types.is_bool_dtype(ann_df[col]):
                ann_df[col] = ann_df[col].astype("category")
        top_annotation = pch.HeatmapAnnotation(df=ann_df, axis=1, plot_legend=True)

    return pch.ClusterMapPlotter(
        data=matrix,
        top_annotation=top_annotation,
        z_score=z_score,
        row_cluster=row_cluster,
        col_cluster=col_cluster,
        cmap=cmap,
        show_rownames=show_rownames,
        show_colnames=show_colnames,
        plot=plot,
        **kwargs,
    )
