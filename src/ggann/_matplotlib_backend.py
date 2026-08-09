"""Opt-in direct Matplotlib rendering for the three primary large plots.

The objects in this module remain :class:`plotnine.ggplot` instances.  A newly
constructed object can draw one proven layout without plotnine's scale/guide
training overhead.  Adding any grammar component invalidates that cached layout
and delegates the complete plot to plotnine, preserving normal composition.
"""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
from typing import Any, Literal

import numpy as np
import pandas as pd
from matplotlib import colormaps
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.colors import BoundaryNorm, LinearSegmentedColormap, ListedColormap
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from mizani.breaks import breaks_extended
from mizani.palettes import area_pal, hue_pal
from plotnine import ggplot


@dataclass(frozen=True)
class _RenderSpec:
    kind: Literal["embedding", "dotplot", "matrixplot"]
    data: pd.DataFrame
    x: str
    y: str
    value: str | None = None
    fraction: str | None = None
    categorical: bool = False
    categories: tuple[Any, ...] = ()
    palette: tuple[str, ...] = ()
    point_size: float = 1.5
    alpha: float = 0.9
    low: str = "#d9d9d9"
    high: str = "#2166ac"
    cmap: str = "viridis"
    size_range: tuple[float, float] = (0.5, 8.0)
    value_label: str = ""
    fraction_label: str = "fraction\nexpressing"


class MatplotlibGGPlot(ggplot):
    """A composable ggplot with an active direct-Matplotlib draw plan.

    Parameters
    ----------
    data : pandas.DataFrame, optional
        Plot data used by the ordinary plotnine grammar.
    mapping : plotnine.mapping.aes, optional
        Plotnine aesthetic mapping.

    Returns
    -------
    MatplotlibGGPlot
        A :class:`plotnine.ggplot` subclass.

    Raises
    ------
    TypeError
        If plotnine rejects ``data`` or ``mapping``.

    Notes
    -----
    Users normally receive this class from a high-level helper with
    ``backend="matplotlib"``. Adding a layer, scale, theme, facet, or annotation
    disables the direct plan on the new plot and restores the full plotnine
    renderer. The original object is unchanged.

    Examples
    --------
    >>> p = plot_embedding(adata, color="cell_type", backend="matplotlib")
    >>> isinstance(p, ggplot)
    True
    """

    _ggann_render_spec: _RenderSpec | None

    @classmethod
    def from_plot(cls, plot: ggplot, spec: _RenderSpec) -> MatplotlibGGPlot:
        """Promote a fully constructed ggplot without copying its data."""
        result = cls.__new__(cls)
        result.__dict__.update(plot.__dict__)
        result._ggann_render_spec = spec
        return result

    @property
    def fast_path_active(self) -> bool:
        """Whether the next draw uses the direct Matplotlib plan."""
        return self._ggann_render_spec is not None

    def __deepcopy__(self, memo: dict[Any, Any]) -> MatplotlibGGPlot:
        """Copy grammar state while retaining immutable plot data by reference."""
        result = self.__class__.__new__(self.__class__)
        memo[id(self)] = result
        shallow = {"data", "figure", "gs", "_build_objs", "_ggann_render_spec"}
        for key, item in self.__dict__.items():
            if key in shallow:
                result.__dict__[key] = item
                memo[id(item)] = item
            else:
                result.__dict__[key] = deepcopy(item, memo)
        return result

    def __iadd__(self, other):
        # Composition means the frozen layout is no longer provably equivalent.
        self._ggann_render_spec = None
        return super().__iadd__(other)

    def draw(self, *, show: bool = False) -> Figure:
        """Draw directly while active, otherwise use plotnine's renderer."""
        spec = self._ggann_render_spec
        if spec is None:
            return super().draw(show=show)
        figure = _draw_spec(spec)
        if show:
            figure.show()
        return figure


def _new_figure(
    *,
    left: float,
    right: float,
    bottom: float,
    top: float,
) -> tuple[Figure, Any]:
    figure = Figure(figsize=(6.4, 4.8))
    FigureCanvasAgg(figure)
    axis = figure.add_axes((left, bottom, right - left, top - bottom))
    axis.grid(False)
    for spine in axis.spines.values():
        spine.set_color("black")
        spine.set_linewidth(1.0)
    return figure, axis


def _categorical_levels(series: pd.Series) -> list[Any]:
    if isinstance(series.dtype, pd.CategoricalDtype):
        return list(series.cat.categories)
    return list(pd.unique(series.dropna()))


def categorical_palette(
    series: pd.Series,
    stored: dict[Any, str] | None,
) -> tuple[tuple[Any, ...], tuple[str, ...]]:
    """Resolve one stable categorical palette without training a plotnine scale."""
    categories = _categorical_levels(series)
    if stored is not None:
        colors = [stored[category] for category in categories]
    else:
        colors = [str(color) for color in hue_pal()(len(categories))]
    return tuple(categories), tuple(colors)


def _embedding_colors(
    spec: _RenderSpec,
) -> tuple[np.ndarray, ListedColormap, BoundaryNorm, list[Line2D]]:
    values = spec.data[spec.value]
    categorical = pd.Categorical(values, categories=spec.categories)
    codes = categorical.codes.astype(np.int16, copy=False)
    palette = list(spec.palette)
    if bool(np.any(codes < 0)):
        codes = np.where(codes < 0, len(palette), codes).astype(np.int16, copy=False)
        palette.append("#7F7F7F")
    cmap = ListedColormap(palette, name="ggann_categories")
    norm = BoundaryNorm(np.arange(len(palette) + 1) - 0.5, len(palette))
    handles = [
        Line2D(
            [],
            [],
            marker="o",
            linestyle="",
            markerfacecolor=color,
            markeredgecolor=color,
            markersize=6,
            label=str(category),
        )
        for category, color in zip(spec.categories, spec.palette)
    ]
    if values.isna().any():
        handles.append(
            Line2D(
                [],
                [],
                marker="o",
                linestyle="",
                markerfacecolor="#7F7F7F",
                markeredgecolor="#7F7F7F",
                markersize=6,
                label="NA",
            )
        )
    return codes, cmap, norm, handles


def _draw_embedding(spec: _RenderSpec) -> Figure:
    right = 0.72 if spec.categorical else (0.84 if spec.value is not None else 0.95)
    figure, axis = _new_figure(left=0.08, right=right, bottom=0.13, top=0.97)
    x = spec.data[spec.x].to_numpy(copy=False)
    y = spec.data[spec.y].to_numpy(copy=False)
    area = np.pi * (spec.point_size + 0.5) ** 2
    if spec.value is None:
        axis.scatter(x, y, s=area, alpha=spec.alpha, linewidths=0, color="#333333")
    elif spec.categorical:
        codes, cmap, norm, handles = _embedding_colors(spec)
        axis.scatter(
            x,
            y,
            s=area,
            alpha=spec.alpha,
            linewidths=0,
            c=codes,
            cmap=cmap,
            norm=norm,
        )
        axis.legend(
            handles=handles,
            title=spec.value,
            loc="center left",
            bbox_to_anchor=(0.75, 0.52),
            bbox_transform=figure.transFigure,
            borderaxespad=0,
            frameon=False,
            fontsize=8,
            title_fontsize=9,
        )
    else:
        cmap = LinearSegmentedColormap.from_list(
            "ggann_expression", [spec.low, spec.high]
        ).with_extremes(bad="#7F7F7F")
        values = spec.data[spec.value].to_numpy(dtype=float, copy=False)
        collection = axis.scatter(
            x,
            y,
            s=area,
            alpha=spec.alpha,
            linewidths=0,
            c=values,
            cmap=cmap,
        )
        colorbar_axis = figure.add_axes((0.88, 0.32, 0.035, 0.48))
        colorbar = figure.colorbar(collection, cax=colorbar_axis)
        _configure_colorbar(colorbar, spec.value, values)
    axis.set_xlabel(spec.x)
    axis.set_ylabel(spec.y)
    axis.set_xticks([])
    axis.set_yticks([])
    return figure


def _axis_categories(series: pd.Series) -> list[Any]:
    levels = _categorical_levels(series)
    observed = set(series.dropna().astype(object))
    return [level for level in levels if level in observed]


def _numeric_limits(values: np.ndarray) -> tuple[float, float]:
    finite = values[np.isfinite(values)]
    if not len(finite):
        return 0.0, 1.0
    low, high = float(finite.min()), float(finite.max())
    return (low, high) if high > low else (low - 0.5, high + 0.5)


def _scale_breaks(values: np.ndarray) -> np.ndarray:
    """Return the same extended continuous breaks used by plotnine scales."""
    finite = values[np.isfinite(values)]
    if not len(finite):
        return np.empty(0, dtype=float)
    low, high = float(finite.min()), float(finite.max())
    breaks = np.asarray(breaks_extended()((low, high)), dtype=float)
    tolerance = np.finfo(float).eps * max(1.0, abs(low), abs(high)) * 8
    return breaks[
        np.isfinite(breaks)
        & (breaks >= low - tolerance)
        & (breaks <= high + tolerance)
    ]


def _dot_sizes(
    values: np.ndarray,
    size_range: tuple[float, float],
    *,
    limits: tuple[float, float] | None = None,
) -> np.ndarray:
    low, high = _numeric_limits(values) if limits is None else limits
    scaled = np.zeros_like(values, dtype=float) if high == low else (values - low) / (high - low)
    sizes = area_pal(size_range)(np.clip(scaled, 0, 1))
    return np.pi * (sizes + 0.5) ** 2


def _configure_colorbar(colorbar, label: str, values: np.ndarray) -> None:
    if colorbar.solids is not None:
        colorbar.solids.set_rasterized(False)
    ticks = _scale_breaks(values)
    if len(ticks):
        colorbar.set_ticks(ticks)
    colorbar.ax.set_title(label, fontsize=9, loc="left", pad=6)


def _set_discrete_axes(axis, x_levels: list[Any], y_levels: list[Any]) -> None:
    axis.set_xticks(np.arange(len(x_levels)), [str(value) for value in x_levels], rotation=45)
    axis.set_yticks(np.arange(len(y_levels)), [str(value) for value in y_levels])
    axis.tick_params(axis="both", labelsize=9)
    axis.set_xlim(-0.5, len(x_levels) - 0.5)
    axis.set_ylim(-0.5, len(y_levels) - 0.5)
    axis.set_xlabel("")
    axis.set_ylabel("")


def _discrete_left_margin(levels: list[Any]) -> float:
    longest = max((len(str(value)) for value in levels), default=0)
    return min(0.55, max(0.16, 0.05 + 0.0105 * longest))


def _draw_dotplot(spec: _RenderSpec) -> Figure:
    x_levels = _axis_categories(spec.data[spec.x])
    y_levels = _axis_categories(spec.data[spec.y])
    figure, axis = _new_figure(
        left=_discrete_left_margin(y_levels),
        right=0.82,
        bottom=0.14,
        top=0.97,
    )
    x_lookup = {value: index for index, value in enumerate(x_levels)}
    y_lookup = {value: index for index, value in enumerate(y_levels)}
    x = np.fromiter((x_lookup[value] for value in spec.data[spec.x]), dtype=float)
    y = np.fromiter((y_lookup[value] for value in spec.data[spec.y]), dtype=float)
    fractions = spec.data[spec.fraction].to_numpy(dtype=float, copy=False)
    values = spec.data[spec.value].to_numpy(dtype=float, copy=False)
    cmap = colormaps[spec.cmap].with_extremes(bad="#7F7F7F")
    collection = axis.scatter(
        x,
        y,
        s=_dot_sizes(fractions, spec.size_range),
        c=values,
        cmap=cmap,
        linewidths=0,
    )
    colorbar_axis = figure.add_axes((0.85, 0.58, 0.035, 0.29))
    colorbar = figure.colorbar(collection, cax=colorbar_axis)
    _configure_colorbar(colorbar, spec.value_label, values)
    fraction_ticks = _scale_breaks(fractions)
    proxy_sizes = _dot_sizes(
        fraction_ticks,
        spec.size_range,
        limits=_numeric_limits(fractions),
    )
    handles = [
        Line2D(
            [],
            [],
            marker="o",
            linestyle="",
            color="#555555",
            markersize=float(np.sqrt(size)),
            label=f"{value:.0%}",
        )
        for value, size in zip(fraction_ticks, proxy_sizes)
    ]
    axis.legend(
        handles=handles,
        title=spec.fraction_label,
        loc="upper left",
        bbox_to_anchor=(0.84, 0.45),
        bbox_transform=figure.transFigure,
        borderaxespad=0,
        frameon=False,
        fontsize=8,
        title_fontsize=9,
    )
    _set_discrete_axes(axis, x_levels, y_levels)
    return figure


def _draw_matrixplot(spec: _RenderSpec) -> Figure:
    x_levels = _axis_categories(spec.data[spec.x])
    y_levels = _axis_categories(spec.data[spec.y])
    figure, axis = _new_figure(
        left=_discrete_left_margin(y_levels),
        right=0.82,
        bottom=0.14,
        top=0.97,
    )
    matrix = (
        spec.data.pivot(index=spec.y, columns=spec.x, values=spec.value)
        .reindex(index=y_levels, columns=x_levels)
        .to_numpy(dtype=float)
    )
    cmap = colormaps[spec.cmap].with_extremes(bad="#7F7F7F")
    mesh = axis.pcolormesh(
        np.arange(len(x_levels) + 1) - 0.5,
        np.arange(len(y_levels) + 1) - 0.5,
        matrix,
        cmap=cmap,
        shading="flat",
        rasterized=False,
    )
    colorbar_axis = figure.add_axes((0.85, 0.39, 0.035, 0.43))
    colorbar = figure.colorbar(mesh, cax=colorbar_axis)
    _configure_colorbar(colorbar, spec.value_label, matrix.ravel())
    _set_discrete_axes(axis, x_levels, y_levels)
    return figure


def _draw_spec(spec: _RenderSpec) -> Figure:
    if spec.kind == "embedding":
        return _draw_embedding(spec)
    if spec.kind == "dotplot":
        return _draw_dotplot(spec)
    return _draw_matrixplot(spec)


def promote_matplotlib_plot(plot: ggplot, **kwargs: Any) -> MatplotlibGGPlot:
    """Attach an immutable direct-render plan to an already complete ggplot."""
    return MatplotlibGGPlot.from_plot(plot, _RenderSpec(**kwargs))


__all__ = ["MatplotlibGGPlot"]
