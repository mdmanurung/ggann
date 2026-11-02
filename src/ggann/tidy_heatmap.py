"""A tidyHeatmap-style verb interface to PyComplexHeatmap.

This is a port of the R package tidyHeatmap (Mangiola and Papenfuss, JOSS 2020).
``tidy_heatmap`` takes a long table (``row``, ``column``, ``value``) or an
``AnnData``; every verb (``annotation_*``, ``layer_*``, ``split_*``) returns a new,
immutable :class:`TidyHeatmap`. Nothing is drawn until ``draw``,
``as_ComplexHeatmap``, ``save`` or Jupyter display.

Like ``plot_clustermap`` this is a grid-backend escape hatch, not a plotnine object.
``PyComplexHeatmap`` is an optional dependency imported lazily, so ``import ggann``
never needs it.
"""

from __future__ import annotations

import contextlib
import dataclasses
import functools
import io
import warnings
from collections.abc import Callable, Mapping, Sequence
from numbers import Real
from typing import Any, NamedTuple

import numpy as np
import pandas as pd
from matplotlib.colors import Colormap, LinearSegmentedColormap, to_hex

from ._aggregate import group_means
from ._compat import SKIP_GGANN_FRAMES, accepts_groupby
from ._expression import densify_frame as _densify
from ._expression import ordered_unique, project_expression, resolve_source
from ._grouping import _drop_missing_groups
from .clustermap import _require_pch
from .publication import _UNIT_TO_INCH, _active_style, _width_in_units, save_publication

__all__ = ["HeatmapData", "TidyHeatmap", "TidyHeatmapSet", "tidy_heatmap"]

_PALETTE_VALUE = ("#440154FF", "#21908CFF", "#fefada")
_SCALES = ("none", "row", "column", "both")
# tidyHeatmap's InputHeatmap prototype: the k-th annotation without a palette takes
# the k-th entry (ColorBrewer qualitative sets; continuous maps by matplotlib name).
_DISCRETE_DEFAULTS = ("Set1", "Set2", "Set3", "Dark2", "Accent", "Pastel2")
_CONTINUOUS_DEFAULTS = ("Spectral_r", "viridis", "magma", "PRGn", "BrBG")
# pch (R pch codes) -> Matplotlib markers; star is pch 8, asterisk is pch 42.
_SYMBOLS = {
    "point": "o",
    "square": "s",
    "diamond": "D",
    "arrow_up": "^",
    "arrow_down": "v",
    "star": (8, 2, 0),
    "asterisk": (6, 2, 0),
}
_RESERVED = {
    "data": "pass the data as the first argument",
    "cmap": "use palette_value",
    "top_annotation": "use the annotation_* verbs",
    "left_annotation": "use the annotation_* verbs",
    "bottom_annotation": "use the annotation_* verbs",
    "right_annotation": "use the annotation_* verbs",
    "row_split": "use split_rows() or group_by()",
    "col_split": "use split_columns() or group_by()",
    "z_score": "use scale=",
    "standard_scale": "use scale=",
}


class HeatmapData(NamedTuple):
    """Matrix in plotted order and the row and column order of the drawn plotter."""

    matrix: pd.DataFrame
    row_order: list
    column_order: list


@dataclasses.dataclass(frozen=True, eq=False)
class _Annot:
    column: str
    kind: str
    axis: str
    palette: Any
    size: float | None
    kwargs: Mapping[str, Any]
    labels_format: Any = None


@dataclasses.dataclass(frozen=True, eq=False)
class _Layer:
    kind: str  # "symbol" or "text"
    symbol: str | None
    cells: pd.DataFrame  # columns: row, column, size, color, text (labels, str)


@dataclasses.dataclass(frozen=True, eq=False)
class _Spec:
    long: pd.DataFrame
    row: str
    column: str
    value: str
    matrix: pd.DataFrame  # transformed, scaled; str labels
    row_attrs: pd.DataFrame
    col_attrs: pd.DataFrame
    unrouted: tuple[str, ...]
    stored_colors: Mapping[str, Mapping]
    palette_value: Any
    kwargs: Mapping[str, Any]
    row_group: str | None = None
    col_group: str | None = None
    palette_grouping: tuple = ()
    group_label_fontsize: float = 8
    show_group_name: bool = True
    group_strip_height: float = 9
    annotations: tuple[_Annot, ...] = ()
    layers: tuple[_Layer, ...] = ()
    row_split: int | None = None
    col_split: int | None = None
    default_show_colnames: bool = True


# --------------------------------------------------------------------------- build


def _entity_order(series: pd.Series, *, sort: bool) -> list:
    """Order of row or column entities: sorted (tidyHeatmap) or as they appear."""
    if isinstance(series.dtype, pd.CategoricalDtype):
        present = set(series.dropna().unique())
        return [c for c in series.cat.categories if c in present]
    unique = list(pd.unique(series.dropna()))
    if not sort:
        return unique
    try:
        return sorted(unique)
    except TypeError:
        return unique


def _scale_robust(m: pd.DataFrame, axis: int) -> pd.DataFrame:
    """tidyHeatmap ``scale_robust``: z-score with ``sd`` (n - 1); sd of 0 or NA -> 1.

    ``axis=1`` scales each row, ``axis=0`` each column.
    """
    mean = m.mean(axis=axis)
    sd = m.std(axis=axis, ddof=1)
    sd = sd.where(sd.notna() & (sd != 0), 1.0)
    other = 1 - axis
    return m.sub(mean, axis=other).div(sd, axis=other)


def _constant_columns(df: pd.DataFrame, key: str, candidates: Sequence[str]) -> list[str]:
    """Candidates that take one value (missing counts as a value) per ``key`` entity."""
    keep = []
    for name in candidates:
        try:
            counts = df.groupby(key, observed=True, sort=False)[name].nunique(dropna=False)
        except TypeError:  # unhashable cells (list columns) cannot be routed
            continue
        if int(counts.max()) <= 1:
            keep.append(name)
    return keep


def _entity_attrs(df, key, columns, order) -> pd.DataFrame:
    """One row per entity (``str`` labelled), with the key itself as a column."""
    first = df.groupby(key, observed=True, sort=False)[list(columns)].first()
    first.index = first.index.astype(str)
    out = first.reindex([str(x) for x in order])
    out.insert(0, key, out.index)
    return out


def _value_cmap(palette_value) -> Colormap:
    if isinstance(palette_value, Colormap):
        return palette_value
    if isinstance(palette_value, str) or len(palette_value) != 3:
        raise ValueError(
            "palette_value must be 3 colours (for example "
            f"{list(_PALETTE_VALUE)}) or a Matplotlib Colormap."
        )
    return LinearSegmentedColormap.from_list("tidy_heatmap", list(palette_value))


def _prepare(
    long: pd.DataFrame,
    row: str,
    column: str,
    value: str,
    *,
    transform,
    scale,
    na,
    sort: bool,
):
    """Validate the long table and return the plotted matrix with its label orders."""
    if transform is not None and not callable(transform):
        raise TypeError("transform must be a callable (for example numpy.log1p) or None.")
    if scale not in _SCALES:
        raise ValueError(f"scale must be one of {list(_SCALES)}, got {scale!r}.")
    if not (na == "raise" or na == "drop" or isinstance(na, (int, float))):
        raise ValueError("na must be 'raise', 'drop' or a number used to fill missing cells.")
    absent = [c for c in (row, column, value) if c not in long.columns]
    if absent:
        raise KeyError(f"Column(s) {absent} not found in the data frame.")
    if len({row, column, value}) != 3:
        raise ValueError("row, column and value must be three different columns.")
    for key in (row, column):
        if int(long[key].isna().sum()):
            raise ValueError(f"Key column {key!r} has missing values.")
    if not pd.api.types.is_numeric_dtype(long[value]) or pd.api.types.is_bool_dtype(long[value]):
        raise ValueError(f"value column {value!r} must be numeric, got {long[value].dtype}.")

    keys = pd.DataFrame({"r": long[row].astype(str), "c": long[column].astype(str)})
    dup = keys.duplicated(keep=False)
    if int(dup.sum()):
        shown = keys[dup].drop_duplicates().head(5).to_numpy().tolist()
        raise ValueError(
            f"Your dataset includes duplicated {row}/{column} pairs ({int(dup.sum())} "
            f"rows, for example {shown}). Remove redundancies before plotting."
        )
    values = long[value].astype(float)
    if transform is not None:
        with np.errstate(all="ignore"):  # NaN and -inf are reported below
            values = pd.Series(np.asarray(transform(values), dtype=float), index=long.index)
        if values.isna().sum() > long[value].isna().sum():
            raise ValueError("You applied a transformation that introduced NaN.")
        if np.isneginf(values).any():
            raise ValueError(
                "You applied a transformation that introduced negative infinite values; "
                "was it log? If so, use log1p."
            )
    if int(values.notna().sum()) == 0:
        raise ValueError("Your dataset does not have any non-missing values.")

    row_order = [str(x) for x in _entity_order(long[row], sort=sort)]
    col_order = [str(x) for x in _entity_order(long[column], sort=sort)]
    wide = (
        pd.DataFrame({"r": keys["r"], "c": keys["c"], "v": values})
        .pivot(index="r", columns="c", values="v")
        .reindex(index=row_order, columns=col_order)
    )
    wide.index.name, wide.columns.name = None, None
    n_missing = int(wide.isna().sum().sum())
    if n_missing:
        if na == "raise":
            raise ValueError(
                f"{n_missing} missing cells in the {row} x {column} matrix (a missing value "
                "or an absent pair). Pass na='drop' to drop those rows, or na=<number> to fill."
            )
        if na == "drop":
            bad = wide.isna().any(axis=1)
            warnings.warn(
                f"dropped {int(bad.sum())} of {len(wide)} {row} entries with missing values",
                UserWarning,
                skip_file_prefixes=SKIP_GGANN_FRAMES,
            )
            wide = wide.loc[~bad]
        else:
            wide = wide.fillna(float(na))
    if wide.empty:
        raise ValueError(f"The matrix to draw is empty ({wide.shape[0]} x {wide.shape[1]}).")
    if scale in ("row", "both"):
        wide = _scale_robust(wide, axis=1)
    if scale in ("column", "both"):
        wide = _scale_robust(wide, axis=0)
    return wide


def _checked_kwargs(kwargs: Mapping[str, Any]) -> dict:
    for name in kwargs:
        if name in _RESERVED:
            raise ValueError(f"{name!r} is managed by tidy_heatmap: {_RESERVED[name]}.")
    return dict(kwargs)


def _group_columns(groupby) -> list[str]:
    keys = groupby.keys
    keys = [keys] if isinstance(keys, str) else list(keys)
    if not all(isinstance(k, str) for k in keys):
        raise ValueError("group_by() keys must be column names.")
    return keys


@accepts_groupby
def tidy_heatmap(
    data,
    row: str | None = None,
    column: str | None = None,
    value: str | None = None,
    *,
    genes: Sequence[str] | None = None,
    group_by: str | None = None,
    layer: str | None = None,
    use_raw: bool | None = None,
    transform: Callable | None = None,
    scale: str = "none",
    palette_value: Sequence[str] | Colormap = _PALETTE_VALUE,
    palette_grouping: Sequence[Sequence[str]] | None = None,
    na: str | float = "raise",
    **kwargs,
) -> TidyHeatmap:
    """Start a tidyHeatmap-style heatmap from a long table or an ``AnnData``.

    A Python port of tidyHeatmap's ``heatmap()``. Every verb on the result returns a
    new :class:`TidyHeatmap`; nothing is drawn until ``draw()``, ``as_ComplexHeatmap()``,
    ``save()`` or Jupyter display.

    * **Long table.** ``data`` is a ``pandas.DataFrame`` with one row per
      (``row``, ``column``) pair. A ``DataFrameGroupBy`` (``df.groupby("lineage")``)
      plays the role of R's ``group_by()``: a grouping column that is constant within
      each row entity splits the rows, one constant within each column entity splits the
      columns, and each split gets a group strip. Rows and columns are sorted (factor
      level order for categoricals) before clustering, as in tidyHeatmap.
    * **AnnData.** Pass ``genes``. Rows are genes and columns are cells, or group means
      when ``group_by`` is given (mean expression, via ``group_means``). ``obs`` columns
      can then be column annotations and ``var`` columns row annotations; with
      ``group_by``, an ``obs`` column is usable when constant within each group, and
      ``n_cells`` is added. Rows and columns keep their given order.

    The matrix is built in this order: ``transform`` on ``value``, the ``na`` policy,
    then ``scale``. ``scale`` uses tidyHeatmap's ``scale_robust``: a z-score with the
    sample standard deviation (n - 1), where a standard deviation of 0 gives a column
    of zeros instead of NaN. ``"both"`` scales rows, then columns.

    Clustering defaults follow ComplexHeatmap (``complete`` linkage, Euclidean
    distance), not PyComplexHeatmap's ``ward``. Value colours interpolate
    ``palette_value`` linearly in RGB between the minimum and maximum of the matrix
    (tidyHeatmap interpolates in LAB space).

    Parameters
    ----------
    data : pandas.DataFrame, pandas.core.groupby.DataFrameGroupBy or anndata.AnnData
        Long table (optionally grouped) or annotated data matrix.
    row, column, value : str
        Long-table column names (tidyHeatmap's ``.row``, ``.column``, ``.value``).
        Not used with ``AnnData``.
    genes : sequence of str, optional
        Genes shown as rows. Required for, and only valid with, ``AnnData``.
    group_by : str, optional
        ``AnnData`` only: obs column whose group means become the columns.
    layer, use_raw : optional
        Mutually exclusive expression source for ``AnnData``.
    transform : callable, optional
        Applied to ``value`` before scaling, for example ``numpy.log1p``. A transform
        that introduces NaN or negative infinity raises.
    scale : {"none", "row", "column", "both"}
        Z-score the rows, the columns or both.
    palette_value : sequence of 3 colours or matplotlib.colors.Colormap
        Heatmap colours, low to high. Default ``#440154FF``, ``#21908CFF``, ``#fefada``.
    palette_grouping : sequence of colour sequences, optional
        Group-strip colours, one sequence per grouping, row grouping first. Unset
        strips are white, as in tidyHeatmap.
    na : {"raise", "drop"} or float
        What to do with a missing value or an absent pair. ``"raise"`` (default) raises,
        ``"drop"`` drops the affected rows with a warning, a number fills them.
    **kwargs
        Passed to ``PyComplexHeatmap.ClusterMapPlotter``, for example
        ``row_cluster=False`` or ``show_rownames=False``. Arguments that the verbs own
        (``cmap``, annotations, splits, scaling) raise ``ValueError``.

    Returns
    -------
    TidyHeatmap
        Immutable builder; chain verbs, then ``draw()`` or ``save()``.

    Raises
    ------
    ImportError
        When a method that draws is called without the ``heatmap`` extra.
    KeyError
        If a named column, gene or layer is missing.
    ValueError
        For duplicated (row, column) pairs, missing values with ``na="raise"``, an
        invalid ``scale`` or ``palette_value``, or incompatible arguments.
    TypeError
        If ``transform`` is not callable.

    Notes
    -----
    Duplicated pairs raise even when their values agree (tidyHeatmap tolerates exact
    duplicates). The input is copied and never modified. Entity labels are compared
    as strings.

    Examples
    --------
    >>> df = pd.DataFrame(
    ...     {"gene": list("abab"), "sample": list("xxyy"), "expr": [1.0, 2.0, 3.0, 5.0]}
    ... )
    >>> hm = tidy_heatmap(df, row="gene", column="sample", value="expr", scale="row")
    >>> hm.get_heatmap_data().matrix.shape  # doctest: +SKIP
    (2, 2)
    """
    kwargs = _checked_kwargs(kwargs)
    if hasattr(data, "obs") and hasattr(data, "var"):
        if any(x is not None for x in (row, column, value)):
            raise ValueError("row, column and value are for long tables; use genes= with AnnData.")
        if genes is None:
            raise ValueError("Pass genes=[...] to build a heatmap from an AnnData.")
        parts = _from_anndata(data, genes, group_by, layer, use_raw)
    else:
        if genes is not None or group_by is not None or layer is not None or use_raw is not None:
            raise ValueError("genes, group_by, layer and use_raw are only for AnnData input.")
        groups: list[str] = []
        if hasattr(data, "obj") and hasattr(data, "keys"):
            groups = _group_columns(data)
            data = data.obj
        if not isinstance(data, pd.DataFrame):
            raise TypeError(
                f"data must be a DataFrame, a DataFrameGroupBy or an AnnData, got {type(data)}."
            )
        if row is None or column is None or value is None:
            raise ValueError("row, column and value are required for a long table.")
        parts = _from_long(data, row, column, value, transform, scale, na)
        parts["groups"] = groups
    if "matrix" not in parts:
        parts["matrix"] = _prepare(
            parts["long"],
            parts["row"],
            parts["column"],
            parts["value"],
            transform=transform,
            scale=scale,
            na=na,
            sort=False,
        )
    groups = parts.pop("groups", [])
    _value_cmap(palette_value)  # fail early on a bad palette
    spec = _Spec(
        palette_value=palette_value,
        kwargs=kwargs,
        palette_grouping=tuple(tuple(p) for p in (palette_grouping or ())),
        **parts,
    )
    if groups:
        spec = _apply_groups(spec, groups)
    return TidyHeatmap(spec)


def _from_long(df, row, column, value, transform, scale, na) -> dict:
    matrix = _prepare(df, row, column, value, transform=transform, scale=scale, na=na, sort=True)
    long = df.copy()
    # route every other column: constant within each row or column entity
    others = [c for c in long.columns if c not in (row, column, value)]
    r_keep = _constant_columns(long, row, others)
    c_keep = _constant_columns(long, column, others)
    row_order = [str(x) for x in _entity_order(long[row], sort=True)]
    col_order = [str(x) for x in _entity_order(long[column], sort=True)]
    row_attrs = _entity_attrs(long, row, r_keep, row_order).loc[matrix.index]
    col_attrs = _entity_attrs(long, column, c_keep, col_order).loc[matrix.columns]
    return {
        "long": long,
        "row": row,
        "column": column,
        "value": value,
        "matrix": matrix,
        "row_attrs": row_attrs,
        "col_attrs": col_attrs,
        "unrouted": tuple(c for c in others if c not in r_keep and c not in c_keep),
        "stored_colors": {},
    }


def _from_anndata(adata, genes, group_by, layer, use_raw) -> dict:
    adata = _drop_missing_groups(adata, group_by)
    genes = ordered_unique(genes)
    kind, lyr = resolve_source(adata, layer, use_raw)
    obs = adata.obs.copy()
    var = (adata.raw.var if kind == "raw" else adata.var).copy()
    if group_by is None:
        projected, genes = project_expression(adata, genes, kind=kind, layer=lyr)
        wide = _densify(projected.ap.to_df(x=genes))[genes].T  # genes x cells
        wide.columns = obs.index.copy()
        column, value = "cell", "expression"
        if not obs.index.is_unique:
            raise ValueError("Cell names (obs_names) must be unique to be heatmap columns.")
        col_attrs = obs
        unrouted: list[str] = []
    else:
        means = group_means(adata, genes, group_by, layer=layer, use_raw=use_raw)
        wide = means.T.loc[genes]  # genes x groups
        column, value = group_by, "mean_expression"
        grouped = obs.groupby(group_by, observed=True, sort=False)
        keep = _constant_columns(obs, group_by, [c for c in obs.columns if c != group_by])
        col_attrs = grouped[keep].first() if keep else pd.DataFrame(index=wide.columns)
        if "n_cells" not in col_attrs:
            col_attrs["n_cells"] = grouped.size()
        unrouted = [c for c in obs.columns if c != group_by and c not in keep]
    row = "gene"
    wide.index.name, wide.columns.name = None, None
    long = wide.reset_index(names=row).melt(id_vars=row, var_name=column, value_name=value)
    col_attrs = col_attrs.copy()
    col_attrs.index = col_attrs.index.astype(str)
    col_attrs = col_attrs.reindex(wide.columns.astype(str))
    col_attrs.insert(0, column, col_attrs.index)
    row_attrs = var.reindex(genes)
    row_attrs.index = row_attrs.index.astype(str)
    row_attrs.insert(0, row, row_attrs.index)
    from ._palette import annotation_colors

    stored = {}
    for name in col_attrs.columns:
        if name in obs and not pd.api.types.is_numeric_dtype(obs[name]):
            colors = annotation_colors(adata, obs[[name]])
            if colors:
                stored[name] = colors[name]
    return {
        "long": long,
        "row": row,
        "column": column,
        "value": value,
        "row_attrs": row_attrs,
        "col_attrs": col_attrs,
        "unrouted": tuple(unrouted),
        "stored_colors": stored,
        "default_show_colnames": group_by is not None,
    }


def _apply_groups(spec: _Spec, groups: Sequence[str]) -> _Spec:
    """Resolve grouping columns to the row or column axis, as tidyHeatmap does."""
    found: dict[str, list[str]] = {"row": [], "column": []}
    for name in groups:
        in_row, in_col = name in spec.row_attrs, name in spec.col_attrs
        if in_row and in_col:
            raise ValueError(
                f"The grouping {name!r} is not specific to rows or columns "
                "(it is constant over the whole table). Maybe you just have one grouping."
            )
        if not in_row and not in_col:
            raise ValueError(
                f"The grouping {name!r} is constant within neither each {spec.row!r} nor "
                f"each {spec.column!r} entity, so it cannot split rows or columns."
            )
        found["row" if in_row else "column"].append(name)
    if max(len(v) for v in found.values()) > 1:
        raise ValueError(
            "At the moment just one grouping per dimension (max 1 row and 1 column) is supported."
        )
    return dataclasses.replace(
        spec,
        row_group=(found["row"] or [spec.row_group])[0],
        col_group=(found["column"] or [spec.col_group])[0],
    )


# -------------------------------------------------------------------------- colours


def _ramp(colors: Sequence, n: int) -> list[str]:
    """``colorRampPalette(colors[1:min(n, len)])(n)``: first ``n`` colours, else interpolate."""
    base = [to_hex(c, keep_alpha=False) for c in list(colors)[:n]]
    if len(base) == n or n == 0:
        return base
    if len(base) == 1:
        return base * n
    cmap = LinearSegmentedColormap.from_list("ramp", base)
    return [to_hex(cmap(x)) for x in np.linspace(0, 1, n)]


def _fill_na(series: pd.Series, what: str) -> pd.Series:
    if not series.isna().any():
        return series
    warnings.warn(
        f"You have NA/NaN values in {what}. These will be replaced with 'NA'.",
        UserWarning,
        skip_file_prefixes=SKIP_GGANN_FRAMES,
    )
    return series.astype(object).where(series.notna(), "NA")


def _levels(series: pd.Series) -> list:
    """Slice or legend order: factor levels if categorical, otherwise sorted."""
    return _entity_order(series, sort=True)


def _discrete_colors(series: pd.Series, palette) -> dict:
    if isinstance(palette, Mapping):
        return {k: palette[k] for k in _levels(series)}
    palette = list(palette)
    present = _levels(series)
    if isinstance(series.dtype, pd.CategoricalDtype):  # tidyHeatmap: palette[seq_len(levels)]
        by_level = dict(zip(series.cat.categories, _ramp(palette, len(series.cat.categories))))
        return {k: by_level[k] for k in present}
    # tidyHeatmap: colorRampPalette(palette)(n_unique) for character and logical columns
    n = len(present)
    if n <= 1:
        return dict(zip(present, palette[:n]))
    cmap = LinearSegmentedColormap.from_list("ramp", [to_hex(c) for c in palette])
    return {k: to_hex(cmap(x)) for k, x in zip(present, np.linspace(0, 1, n))}


def _continuous_cmap(palette) -> Colormap:
    import matplotlib

    if isinstance(palette, Colormap):
        return palette
    if isinstance(palette, str):
        return matplotlib.colormaps[palette]
    return LinearSegmentedColormap.from_list("annotation", [to_hex(c) for c in palette])


# ------------------------------------------------------------------------- drawing


@functools.cache
def _numeric_annotation():
    """Bar annotation that prints each value, the closest PyComplexHeatmap has to anno_numeric."""
    pch = _require_pch()

    class anno_numeric(pch.anno_barplot):
        def __init__(self, df=None, labels_format=None, **kwargs):
            self._format = labels_format
            super().__init__(df, **kwargs)

        def plot(self, ax=None, axis=1):
            ax = super().plot(ax=ax, axis=axis)
            fmt = self._format
            label = fmt if callable(fmt) else (lambda x: (fmt or "{:.1f}").format(x))
            extent = ax.get_window_extent()
            cell = (extent.height if axis == 0 else extent.width) / max(self.nrows, 1)
            size = min(8.0, 0.75 * cell * 72 / ax.figure.dpi)
            for pos, v in zip(np.arange(0.5, self.nrows, 1), self.plot_data.iloc[:, 0]):
                if axis == 1:
                    ax.text(pos, self.ylim[0], label(v), rotation=90, ha="center", va="bottom")
                else:
                    ax.text(self.ylim[0], pos, label(v), ha="left", va="center")
            for text in ax.texts:
                text.set_fontsize(size)
                text.set_zorder(20)
            return ax

    return anno_numeric


def _annotations_for(spec: _Spec, pch, axis: str):
    """Group strips and ``annotation_*`` verbs for ``axis`` as one HeatmapAnnotation."""
    attrs = spec.row_attrs if axis == "row" else spec.col_attrs
    labels = list(spec.matrix.index if axis == "row" else spec.matrix.columns)
    horizontal = axis == "column"
    items: dict[str, Any] = {}

    group = spec.row_group if axis == "row" else spec.col_group
    if group is not None:
        series = _fill_na(attrs[group].reindex(labels), f"your {axis} grouping column")
        index = 0 if axis == "row" else (1 if spec.row_group else 0)
        palettes = spec.palette_grouping
        palette = palettes[index] if index < len(palettes) else None
        levels = _levels(series)
        if palette is None:
            colors = dict.fromkeys(levels, "#ffffff")
        else:
            colors = dict(zip(levels, _ramp(palette, len(levels))))
        height = spec.group_strip_height * 25.4 / 72  # points -> mm
        name = group if spec.show_group_name else " " * (1 if axis == "row" else 2)
        items[name] = pch.anno_simple(
            pd.Series(pd.Categorical(series, categories=levels), index=series.index, name=name),
            colors=colors,
            add_text=True,
            legend=False,
            height=height,
            text_kws={"fontsize": spec.group_label_fontsize},
        )

    n_discrete = n_continuous = 0
    for ann in spec.annotations:
        if ann.axis != axis:
            continue
        series = attrs[ann.column].reindex(labels)
        series.name = ann.column
        size = ann.size
        extra = dict(ann.kwargs)
        numeric = pd.api.types.is_numeric_dtype(series) and not pd.api.types.is_bool_dtype(series)
        if ann.kind == "tile" and not numeric:
            series = _fill_na(series, f"annotation {ann.column!r}")
            if ann.palette is not None:
                palette = ann.palette
            elif ann.column in spec.stored_colors:
                palette = spec.stored_colors[ann.column]
            else:
                if n_discrete >= len(_DISCRETE_DEFAULTS):
                    raise ValueError(
                        "Your discrete annotation columns outnumber the default palettes "
                        f"({len(_DISCRETE_DEFAULTS)}); pass palette= for some of them."
                    )
                palette = _mpl_colors(_DISCRETE_DEFAULTS[n_discrete])
                n_discrete += 1
            colors = _discrete_colors(series, palette)
            data = pd.Series(
                pd.Categorical(series, categories=list(colors)), index=series.index, name=ann.column
            )
            items[ann.column] = pch.anno_simple(data, colors=colors, height=size, **extra)
        elif ann.kind == "tile":
            if ann.palette is not None:
                cmap = _continuous_cmap(ann.palette)
            else:
                if n_continuous >= len(_CONTINUOUS_DEFAULTS):
                    raise ValueError(
                        "Your continuous annotation columns outnumber the default palettes "
                        f"({len(_CONTINUOUS_DEFAULTS)}); pass palette= for some of them."
                    )
                cmap = _continuous_cmap(_CONTINUOUS_DEFAULTS[n_continuous])
                n_continuous += 1
            items[ann.column] = pch.anno_simple(
                series.astype(float), cmap=cmap, height=size, **extra
            )
        else:
            values = series.astype(float)
            color = ann.palette if isinstance(ann.palette, str) else None
            if ann.kind == "point":
                anno = pch.anno_scatterplot(values, colors=color or "black", height=size, **extra)
            elif ann.kind == "line":
                anno = pch.anno_lineplot(values, colors=color or "black", height=size, **extra)
            else:
                lo, hi = min(0.0, float(values.min())), max(0.0, float(values.max()))
                pad = 0.05 * ((hi - lo) or 1.0)
                extra.setdefault("ylim", (lo - pad if lo < 0 else 0.0, hi + pad if hi > 0 else 0.0))
                if ann.kind == "numeric":
                    extra["labels_format"] = ann.labels_format
                    cls, fill = _numeric_annotation(), "#b3b3b3"
                else:
                    cls, fill = pch.anno_barplot, "#cccccc"
                extra.setdefault("colors", color or fill)
                if len(values) > 100:  # outlines would swamp thin bars
                    extra.setdefault("edgecolor", "none")
                anno = cls(values, height=size, **extra)
            items[ann.column] = anno
    if not items:
        return None
    with contextlib.redirect_stdout(io.StringIO()):
        return pch.HeatmapAnnotation(
            axis=1 if horizontal else 0, verbose=0, plot_legend=True, **items
        )


def _mpl_colors(name: str) -> list[str]:
    import matplotlib

    return [to_hex(c) for c in matplotlib.colormaps[name].colors]


def _default_figsize(spec: _Spec) -> tuple[float, float]:
    n_rows, n_cols = spec.matrix.shape
    top = sum(a.axis == "column" for a in spec.annotations) + bool(spec.col_group)
    left = sum(a.axis == "row" for a in spec.annotations) + bool(spec.row_group)
    width = min(14.0, max(5.0, 3.5 + 0.3 * n_cols + 0.5 * left)) + 2.2  # + legends
    height = min(14.0, max(4.0, 2.5 + 0.28 * n_rows + 0.5 * top))
    return width, height


def _cut(matrix: pd.DataFrame, k: int, method: str, metric: str):
    """``cutree(hclust(dist(matrix)), k)``: slice labels 1..k in dendrogram leaf order.

    PyComplexHeatmap's own ``row_split=<int>`` fails on a singleton slice, so the cut is
    done here and passed as an explicit split with an explicit slice order.
    """
    from scipy.cluster.hierarchy import fcluster, leaves_list, linkage

    tree = linkage(matrix.to_numpy(), method=method, metric=metric)
    cut = fcluster(tree, t=k, criterion="maxclust")
    by_leaf = pd.Series(cut[leaves_list(tree)]).drop_duplicates().tolist()
    names = {old: new for new, old in enumerate(by_leaf, start=1)}
    return pd.Series([names[c] for c in cut], index=matrix.index, dtype=object), list(
        names.values()
    )


def _plotter(spec: _Spec, *, plot: bool, defaults: Mapping[str, Any] | None = None):
    pch = _require_pch()
    matrix = spec.matrix
    kwargs = dict(spec.kwargs)
    lo, hi = float(np.nanmin(matrix.to_numpy())), float(np.nanmax(matrix.to_numpy()))
    if lo == hi:  # tidyHeatmap issue 40: one value only
        lo, hi = lo - 1, hi + 1
    n_rows, n_cols = matrix.shape
    for flag, axis, n in (("row_cluster", "row", n_rows), ("col_cluster", "column", n_cols)):
        if kwargs.get(flag, True) and n < 2:
            raise ValueError(
                f"{axis} clustering needs at least 2 {axis}s, got {n}; pass {flag}=False."
            )

    args: dict[str, Any] = {
        "cmap": _value_cmap(spec.palette_value),
        "vmin": lo,
        "vmax": hi,
        "label": spec.value,
        "xlabel": spec.column,
        "ylabel": spec.row,
        "row_cluster_method": "complete",
        "col_cluster_method": "complete",
        "show_rownames": True,
        "show_colnames": spec.default_show_colnames,
        "row_dendrogram": kwargs.get("row_cluster", True),
        "col_dendrogram": kwargs.get("col_cluster", True),
        "yticklabels_kws": {"labelsize": min(12.0, 320 / n_rows)},
        "xticklabels_kws": {"labelsize": min(12.0, 320 / n_cols)},
        "verbose": 0,
    }
    if spec.row_group is not None:
        series = _fill_na(spec.row_attrs[spec.row_group].reindex(matrix.index), "row groups")
        args["row_split"] = series.astype(object)
        args["row_split_order"] = _levels(series)
    elif spec.row_split:
        args["row_split"], args["row_split_order"] = _cut(
            matrix,
            spec.row_split,
            kwargs.get("row_cluster_method", "complete"),
            kwargs.get("row_cluster_metric", "euclidean"),
        )
    if spec.col_group is not None:
        series = _fill_na(spec.col_attrs[spec.col_group].reindex(matrix.columns), "column groups")
        args["col_split"] = series.astype(object)
        args["col_split_order"] = _levels(series)
    elif spec.col_split:
        args["col_split"], args["col_split_order"] = _cut(
            matrix.T,
            spec.col_split,
            kwargs.get("col_cluster_method", "complete"),
            kwargs.get("col_cluster_metric", "euclidean"),
        )
    args.update(defaults or {})
    args.update(kwargs)

    with warnings.catch_warnings(), contextlib.redirect_stdout(io.StringIO()):
        warnings.filterwarnings("ignore", message="Clustering large matrix with scipy")
        warnings.filterwarnings("ignore", message=".*set_bad.*", category=PendingDeprecationWarning)
        result = pch.ClusterMapPlotter(
            data=matrix,
            top_annotation=_annotations_for(spec, pch, "column"),
            left_annotation=_annotations_for(spec, pch, "row"),
            plot=plot,
            **args,
        )
    if style := _active_style():
        result._ggann_publication_style = style
    return result


def _draw_layers(spec: _Spec, plotter, row_order) -> None:
    """Draw the ``layer_*`` symbols and text on the heatmap axes of a drawn plotter."""
    if not spec.layers:
        return
    where_r = {lab: (i, k) for i, labs in enumerate(row_order) for k, lab in enumerate(labs)}
    where_c = {
        lab: (j, k) for j, labs in enumerate(plotter.col_order) for k, lab in enumerate(labs)
    }
    for layer in spec.layers:
        cells = layer.cells
        pos_r = cells["row"].map(where_r)
        pos_c = cells["column"].map(where_c)
        ok = pos_r.notna() & pos_c.notna()
        for idx in np.flatnonzero(ok.to_numpy()):
            i, y = pos_r.iloc[idx]
            j, x = pos_c.iloc[idx]
            ax = plotter.heatmap_axes[i, j]
            limits = ax.get_xlim(), ax.get_ylim()
            cell = cells.iloc[idx]
            if layer.kind == "text":
                ax.text(
                    x + 0.5,
                    y + 0.5,
                    cell["text"],
                    ha="center",
                    va="center",
                    fontsize=cell["size"],
                    color=cell["color"],
                    zorder=20,
                )
            else:
                color, marker = cell["color"], _SYMBOLS[layer.symbol]
                # line markers (star, asterisk) take their colour from the face colour only
                colors = (
                    {"facecolors": color}
                    if isinstance(marker, tuple)
                    else {
                        "facecolors": color,
                        "edgecolors": color,
                    }
                )
                ax.scatter(
                    [x + 0.5],
                    [y + 0.5],
                    s=(float(cell["size"]) * 72 / 25.4) ** 2,
                    marker=marker,
                    linewidths=0.8,
                    zorder=20,
                    **colors,
                )
            ax.set_xlim(limits[0])
            ax.set_ylim(limits[1])


def _export_figsize(width, height, units) -> tuple[float, float]:
    """Canvas size in inches, validated like ``save_publication`` before anything is drawn."""
    width_in, _ = _width_in_units(width, units)
    if isinstance(height, bool) or not isinstance(height, Real) or not float(height) > 0:
        raise ValueError(f"height must be a positive number, got {height!r}.")
    return width_in, float(height) * _UNIT_TO_INCH[units]


def _figure(figsize):
    import matplotlib.pyplot as plt

    return plt.figure(figsize=figsize)


def _png(draw) -> bytes:
    import matplotlib.pyplot as plt

    fig = None
    try:
        draw()
        fig = plt.gcf()
        buffer = io.BytesIO()
        fig.savefig(buffer, format="png", dpi=150, bbox_inches="tight")
        return buffer.getvalue()
    finally:
        if fig is not None:
            plt.close(fig)


# ---------------------------------------------------------------------------- class


def _as_list(columns) -> list[str]:
    return [columns] if isinstance(columns, str) else list(columns)


class TidyHeatmap:
    """Immutable heatmap builder returned by :func:`tidy_heatmap`.

    Verbs (``annotation_*``, ``layer_*``, ``split_*``, ``annotation_group``) return a new
    ``TidyHeatmap``; the original is never changed. ``draw``, ``as_ComplexHeatmap``,
    ``save``, ``get_heatmap_data`` and Jupyter display render it with PyComplexHeatmap.
    Create one with :func:`tidy_heatmap`.

    Parameters
    ----------
    spec : object
        Internal state. Do not construct directly.

    Returns
    -------
    TidyHeatmap
        The builder instance.

    Raises
    ------
    AttributeError
        On any attempt to set an attribute: instances are immutable.

    Notes
    -----
    ``a + b`` places two heatmaps side by side (see :class:`TidyHeatmapSet`).

    Examples
    --------
    >>> hm = tidy_heatmap(df, row="gene", column="sample", value="expr")  # doctest: +SKIP
    >>> hm.annotation_tile("group").split_rows(2).draw()  # doctest: +SKIP
    """

    __slots__ = ("_spec",)

    def __init__(self, spec: _Spec):
        object.__setattr__(self, "_spec", spec)

    def __setattr__(self, name, value):
        raise AttributeError("TidyHeatmap is immutable; verbs return a new TidyHeatmap.")

    def __delattr__(self, name):
        raise AttributeError("TidyHeatmap is immutable; verbs return a new TidyHeatmap.")

    def __repr__(self) -> str:
        s = self._spec
        notes = [f"{s.matrix.shape[0]} {s.row} x {s.matrix.shape[1]} {s.column}"]
        if s.annotations:
            notes.append(f"{len(s.annotations)} annotation(s)")
        if s.layers:
            notes.append(f"{len(s.layers)} layer(s)")
        return f"TidyHeatmap({', '.join(notes)})"

    def _with(self, **changes) -> TidyHeatmap:
        return TidyHeatmap(dataclasses.replace(self._spec, **changes))

    # ---- annotations -------------------------------------------------------

    def _annotate(self, kind, columns, palette, size, kwargs, labels_format=None) -> TidyHeatmap:
        spec = self._spec
        added = []
        for name in _as_list(columns):
            axes = [
                ax
                for ax, attrs in (("row", spec.row_attrs), ("column", spec.col_attrs))
                if name in attrs.columns
            ]
            if not axes:
                if name in spec.unrouted:
                    raise ValueError(
                        f"Your annotation {name!r} is not unique to vertical nor horizontal "
                        f"dimensions: it is constant within neither each {spec.row!r} nor "
                        f"each {spec.column!r} entity."
                    )
                raise KeyError(f"Annotation column {name!r} not found in the data.")
            for axis in axes:
                if any(
                    a.column == name and a.axis == axis for a in spec.annotations + tuple(added)
                ):
                    raise ValueError(f"{name!r} is already an annotation on the {axis}s.")
                attrs = spec.row_attrs if axis == "row" else spec.col_attrs
                if kind != "tile" and not pd.api.types.is_numeric_dtype(attrs[name]):
                    raise ValueError(
                        f"annotation_{kind} needs a numeric column; {name!r} is {attrs[name].dtype}."
                    )
                added.append(_Annot(name, kind, axis, palette, size, dict(kwargs), labels_format))
        return self._with(annotations=spec.annotations + tuple(added))

    def annotation_tile(self, column, palette=None, size=None, **kwargs) -> TidyHeatmap:
        """Add a colour-tile annotation (tidyHeatmap ``annotation_tile``).

        The column goes to the rows or the columns by routing: constant within each row
        entity (rows) or each column entity (columns). ``palette`` is a colour sequence
        or ``{value: colour}`` for a discrete column, a colour sequence or Matplotlib
        colormap (or its name) for a continuous one. Without it, discrete tiles use the
        tidyHeatmap defaults in order (Set1, Set2, Set3, Dark2, Accent, Pastel2), reuse
        colours stored in ``adata.uns`` for ``AnnData`` input, and continuous tiles use
        Spectral (reversed), viridis, magma, PRGn, BrBG. ``size`` is the strip height
        (width for row annotations) in millimetres; ``**kwargs`` go to
        ``PyComplexHeatmap.anno_simple``. A list of column names adds one tile each.
        """
        return self._annotate("tile", column, palette, size, kwargs)

    def annotation_point(self, column, palette=None, size=None, **kwargs) -> TidyHeatmap:
        """Add a point (scatter) annotation for a numeric column (``annotation_point``).

        ``size`` is the annotation height in millimetres; a single colour string as
        ``palette`` sets the point colour (default black). ``**kwargs`` go to
        ``PyComplexHeatmap.anno_scatterplot``.
        """
        return self._annotate("point", column, palette, size, kwargs)

    def annotation_bar(self, column, palette=None, size=None, **kwargs) -> TidyHeatmap:
        """Add a bar annotation for a numeric column (``annotation_bar``).

        Bars start at 0. ``size`` is the annotation height in millimetres; a single
        colour string as ``palette`` sets the fill (default ``#cccccc``). ``**kwargs``
        go to ``PyComplexHeatmap.anno_barplot``.
        """
        return self._annotate("bar", column, palette, size, kwargs)

    def annotation_line(self, column, palette=None, size=None, **kwargs) -> TidyHeatmap:
        """Add a line annotation for a numeric column (``annotation_line``).

        ``size`` is the annotation height in millimetres. ``**kwargs`` go to
        ``PyComplexHeatmap.anno_lineplot``.
        """
        return self._annotate("line", column, palette, size, kwargs)

    def annotation_numeric(
        self, column, palette=None, size=None, labels_format=None, **kwargs
    ) -> TidyHeatmap:
        """Add a bar annotation that prints each value (``annotation_numeric``).

        Partly supported: PyComplexHeatmap has no ``anno_numeric``, so this is a bar
        annotation with the value written inside each cell (grey bars, as in tidyHeatmap).
        ``labels_format`` is a callable or a ``str.format`` template (default
        ``"{:.1f}"``).
        """
        return self._annotate("numeric", column, palette, size, kwargs, labels_format)

    def annotation_group(
        self,
        *columns: str,
        palette_grouping=None,
        group_label_fontsize: float = 8,
        show_group_name: bool = True,
        group_strip_height: float = 9,
    ) -> TidyHeatmap:
        """Split by grouping columns and draw group strips (``annotation_group``).

        The equivalent of calling ``group_by()`` on the input: each column is routed to
        the rows or the columns and splits that axis. ``group_strip_height`` is in
        points. Fully supported except ``show_group_name=False``, which blanks the strip
        label instead of removing it.
        """
        spec = self._spec
        updates: dict[str, Any] = {
            "group_label_fontsize": group_label_fontsize,
            "show_group_name": show_group_name,
            "group_strip_height": group_strip_height,
        }
        if palette_grouping is not None:
            updates["palette_grouping"] = tuple(tuple(p) for p in palette_grouping)
        spec = dataclasses.replace(spec, **updates)
        if columns:
            spec = _apply_groups(spec, list(columns))
        return TidyHeatmap(spec)

    # ---- layers ------------------------------------------------------------

    def _layer(self, kind, symbol, where, size, color, label=None) -> TidyHeatmap:
        spec = self._spec
        frame = spec.long
        for key, attrs in ((spec.row, spec.row_attrs), (spec.column, spec.col_attrs)):
            new = [c for c in attrs.columns if c not in frame.columns]
            if new:
                frame = frame.join(attrs[new], on=frame[key].astype(str))
        if where is None:
            sel = frame
        elif callable(where):
            sel = frame[np.asarray(where(frame), dtype=bool)]
        elif isinstance(where, str):
            sel = frame.query(where)
        else:
            raise TypeError("where must be a pandas query string, a callable or None.")
        n_rows, n_cols = spec.matrix.shape
        if kind == "text":
            default_size = min(12.0, 320 / max(n_rows, n_cols))
        else:
            default_size = min(3.0, 80 / max(n_rows, n_cols))

        def resolve(arg, default, name):
            if arg is None:
                return pd.Series(default, index=sel.index)
            if callable(arg):
                return pd.Series(np.asarray(arg(sel)), index=sel.index)
            if isinstance(arg, str) and arg in sel.columns:
                return sel[arg]
            if name == "size" and isinstance(arg, str):
                raise KeyError(f"size column {arg!r} not found in the long table.")
            return pd.Series(arg, index=sel.index)

        cells = pd.DataFrame(
            {
                "row": sel[spec.row].astype(str),
                "column": sel[spec.column].astype(str),
                "size": resolve(size, default_size, "size"),
                "color": resolve(color, "#000000" if kind == "text" else "#161616", "color"),
                "text": sel[label].astype(str) if label is not None else "",
            }
        )
        return self._with(layers=spec.layers + (_Layer(kind, symbol, cells),))

    def layer_point(self, where=None, *, size=None, color=None) -> TidyHeatmap:
        """Mark cells with filled circles (``layer_point``, R symbol 21).

        ``where`` selects cells from the long table: a pandas ``query`` string (for
        example ``"padj < 0.05"``), a callable returning a boolean mask, or ``None`` for
        every cell. ``size`` (millimetres) and ``color`` are a constant, a column name or a
        callable over the selected rows. Defaults follow tidyHeatmap: size
        ``min(3, 80 / max(n_rows, n_columns))`` mm, colour ``#161616``. The long table
        also exposes annotation columns (``obs`` or ``var`` fields for ``AnnData`` input).
        """
        return self._layer("symbol", "point", where, size, color)

    def layer_square(self, where=None, *, size=None, color=None) -> TidyHeatmap:
        """Mark cells with squares (``layer_square``, R symbol 22); arguments as ``layer_point``."""
        return self._layer("symbol", "square", where, size, color)

    def layer_diamond(self, where=None, *, size=None, color=None) -> TidyHeatmap:
        """Mark cells with diamonds (``layer_diamond``, R symbol 23); arguments as ``layer_point``."""
        return self._layer("symbol", "diamond", where, size, color)

    def layer_arrow_up(self, where=None, *, size=None, color=None) -> TidyHeatmap:
        """Mark cells with up triangles (``layer_arrow_up``, R symbol 24); arguments as ``layer_point``."""
        return self._layer("symbol", "arrow_up", where, size, color)

    def layer_arrow_down(self, where=None, *, size=None, color=None) -> TidyHeatmap:
        """Mark cells with down triangles (``layer_arrow_down``, R symbol 25); arguments as ``layer_point``."""
        return self._layer("symbol", "arrow_down", where, size, color)

    def layer_star(self, where=None, *, size=None, color=None) -> TidyHeatmap:
        """Mark cells with 8-spoke stars (``layer_star``, R symbol 8); arguments as ``layer_point``."""
        return self._layer("symbol", "star", where, size, color)

    def layer_asterisk(self, where=None, *, size=None, color=None) -> TidyHeatmap:
        """Mark cells with 6-spoke asterisks (``layer_asterisk``, R symbol 42); arguments as ``layer_point``."""
        return self._layer("symbol", "asterisk", where, size, color)

    def layer_text(self, where=None, *, label, size=None, color=None) -> TidyHeatmap:
        """Write a label in cells (``layer_text``).

        ``label`` is the long-table column to print. ``where``, ``size`` (here a font size
        in points, default ``min(12, 320 / max(n_rows, n_columns))``) and ``color``
        (default black) work as in ``layer_point``.
        """
        if label not in self._spec.long.columns:
            attrs = list(self._spec.row_attrs.columns) + list(self._spec.col_attrs.columns)
            if label not in attrs:
                raise KeyError(f"Label column {label!r} not found in the long table.")
        return self._layer("text", None, where, size, color, label)

    # ---- splits ------------------------------------------------------------

    def _split(self, axis: str, number_of_groups: int) -> TidyHeatmap:
        spec = self._spec
        n = spec.matrix.shape[0 if axis == "row" else 1]
        if isinstance(number_of_groups, bool) or not isinstance(
            number_of_groups, (int, np.integer)
        ):
            raise TypeError("number_of_groups must be an integer.")
        if not 2 <= number_of_groups <= n:
            raise ValueError(f"number_of_groups must be between 2 and {n}, got {number_of_groups}.")
        if spec.kwargs.get("row_cluster" if axis == "row" else "col_cluster", True) is False:
            raise ValueError(f"split_{axis}s cuts the dendrogram, but clustering is switched off.")
        if (spec.row_group if axis == "row" else spec.col_group) is not None:
            raise ValueError(
                f"The {axis}s are already split by a group_by() column; "
                f"split_{axis}s would replace it."
            )
        return self._with(**{f"{axis[:3]}_split": int(number_of_groups)})

    def split_rows(self, number_of_groups: int) -> TidyHeatmap:
        """Cut the row dendrogram into ``number_of_groups`` slices (``split_rows``).

        Uses the heatmap's clustering method and metric (``complete`` and Euclidean by
        default). Raises if rows are already split by ``group_by()`` or clustering is off.
        """
        return self._split("row", number_of_groups)

    def split_columns(self, number_of_groups: int) -> TidyHeatmap:
        """Cut the column dendrogram into ``number_of_groups`` slices (``split_columns``)."""
        return self._split("column", number_of_groups)

    # ---- drawing -----------------------------------------------------------

    def draw(self, figsize: tuple[float, float] | None = None):
        """Draw into a new Matplotlib figure and return the ``ClusterMapPlotter``.

        ``figsize`` is in inches; PyComplexHeatmap lays out at the size of the figure
        when it is drawn, so pass the export size. The figure stays open, reachable
        as ``plotter.ax.figure``.
        """
        fig = _figure(figsize or _default_figsize(self._spec))
        try:
            plotter = _plotter(self._spec, plot=True)
            _draw_layers(self._spec, plotter, plotter.row_order)
        except BaseException:
            import matplotlib.pyplot as plt

            plt.close(fig)
            raise
        return plotter

    def as_ComplexHeatmap(self, plot: bool = True, figsize: tuple[float, float] | None = None):
        """Return the ``PyComplexHeatmap.ClusterMapPlotter`` (tidyHeatmap ``as_ComplexHeatmap``).

        With ``plot=True`` (default) it is drawn, including the ``layer_*`` marks, as in
        :meth:`draw`. With ``plot=False`` it is returned unplotted and without layers,
        ready for ``composite`` or ``.plot()``.
        """
        if plot:
            return self.draw(figsize)
        return _plotter(self._spec, plot=False)

    def get_heatmap_data(self) -> HeatmapData:
        """Return the matrix in plotted order with the row and column order.

        Draws the heatmap on a temporary figure (closed afterwards) and reads the order
        the plotter used. The matrix is the transformed and scaled matrix. With splits,
        the order runs slice by slice.
        """
        import matplotlib.pyplot as plt

        plotter = self.draw()
        try:
            rows = [r for part in plotter.row_order for r in part]
            cols = [c for part in plotter.col_order for c in part]
            return HeatmapData(plotter.data2d.loc[rows, cols].copy(), rows, cols)
        finally:
            plt.close(plotter.ax.figure)

    def save(
        self,
        path,
        *,
        width="single-column",
        height: float,
        units: str = "mm",
        formats=None,
        dpi: int = 600,
        background: str = "white",
    ):
        """Save with :func:`ggann.save_publication` (same arguments and contract).

        The heatmap is drawn first at the export size, because PyComplexHeatmap lays out
        at the size of the figure it is created on.
        """
        import matplotlib.pyplot as plt

        plotter = self.draw(figsize=_export_figsize(width, height, units))
        try:
            return save_publication(
                plotter,
                path,
                width=width,
                height=height,
                units=units,
                formats=formats,
                dpi=dpi,
                background=background,
            )
        finally:
            plt.close(plotter.ax.figure)

    def _repr_png_(self) -> bytes:
        return _png(self.draw)

    def __add__(self, other) -> TidyHeatmapSet:
        if isinstance(other, TidyHeatmap):
            return TidyHeatmapSet((self, other))
        if isinstance(other, TidyHeatmapSet):
            return TidyHeatmapSet((self, *other.heatmaps))
        return NotImplemented


class TidyHeatmapSet:
    """Heatmaps placed side by side with ``+`` (tidyHeatmap's ``+.InputHeatmap``).

    Uses ``PyComplexHeatmap.composite``. All heatmaps must have the same rows; the
    first one fixes the row order and the others follow it. PyComplexHeatmap does not
    align the heatmap bodies when annotation heights differ, so give the heatmaps the
    same number of column annotations. Legends are collected at the right edge.

    Parameters
    ----------
    heatmaps : sequence of TidyHeatmap
        Heatmaps in left-to-right order.

    Returns
    -------
    TidyHeatmapSet
        Immutable collection with ``draw``, ``save`` and Jupyter display.

    Raises
    ------
    ValueError
        If the heatmaps do not share the same rows.

    Notes
    -----
    Row splits and ``group_by()`` row groups follow the first heatmap only if the others
    use ``row_cluster=False``; otherwise PyComplexHeatmap reuses the first order.

    Examples
    --------
    >>> side_by_side = hm_a + hm_b  # doctest: +SKIP
    >>> fig = side_by_side.draw()  # doctest: +SKIP
    """

    __slots__ = ("heatmaps",)

    def __init__(self, heatmaps: Sequence[TidyHeatmap]):
        heatmaps = tuple(heatmaps)
        first = set(heatmaps[0]._spec.matrix.index)
        if any(set(h._spec.matrix.index) != first for h in heatmaps[1:]):
            raise ValueError("Heatmaps added with + must have the same rows.")
        object.__setattr__(self, "heatmaps", heatmaps)

    def __setattr__(self, name, value):
        raise AttributeError("TidyHeatmapSet is immutable.")

    def __repr__(self) -> str:
        return f"TidyHeatmapSet({len(self.heatmaps)} heatmaps)"

    def __add__(self, other) -> TidyHeatmapSet:
        if isinstance(other, TidyHeatmap):
            return TidyHeatmapSet((*self.heatmaps, other))
        if isinstance(other, TidyHeatmapSet):
            return TidyHeatmapSet((*self.heatmaps, *other.heatmaps))
        return NotImplemented

    def draw(self, figsize: tuple[float, float] | None = None):
        """Draw all heatmaps in one new figure and return the Matplotlib ``Figure``."""
        import matplotlib.pyplot as plt

        pch = _require_pch()
        sizes = [_default_figsize(h._spec) for h in self.heatmaps]
        fig = _figure(figsize or (sum(w for w, _ in sizes) + 1.0, max(h for _, h in sizes)))
        try:
            # only the first heatmap clusters rows; the others follow its order
            follow = {"row_dendrogram": False, "show_rownames": False}
            plotters = [
                _plotter(h._spec, plot=False, defaults=follow if i else None)
                for i, h in enumerate(self.heatmaps)
            ]
            with warnings.catch_warnings(), contextlib.redirect_stdout(io.StringIO()):
                warnings.filterwarnings("ignore", message="Clustering large matrix with scipy")
                pch.composite(plotters, ax=plt.gca(), axis=1, verbose=0)
            for h, plotter in zip(self.heatmaps, plotters):
                _draw_layers(h._spec, plotter, plotters[0].row_order)
        except BaseException:
            plt.close(fig)
            raise
        return fig

    def save(
        self,
        path,
        *,
        width="single-column",
        height: float,
        units: str = "mm",
        formats=None,
        dpi: int = 600,
        background: str = "white",
    ):
        """Save with :func:`ggann.save_publication`; the figure is drawn at the export size."""
        import matplotlib.pyplot as plt

        fig = self.draw(figsize=_export_figsize(width, height, units))
        try:
            return save_publication(
                fig,
                path,
                width=width,
                height=height,
                units=units,
                formats=formats,
                dpi=dpi,
                background=background,
            )
        finally:
            plt.close(fig)

    def _repr_png_(self) -> bytes:
        return _png(self.draw)
