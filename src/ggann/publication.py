"""Journal-oriented styling, colour vocabularies, and exact figure export.

The public API in this module is intentionally independent of any one journal.
It adapts the useful ideas of centralised physical sizing, editable vector
output, and restrained typography from CNSPlots without depending on, or
copying implementation from, that project.
"""

from __future__ import annotations

import colorsys
import math
import warnings
from collections.abc import Iterator, Mapping, Sequence
from contextlib import contextmanager
from contextvars import ContextVar
from copy import deepcopy
from dataclasses import asdict, dataclass, replace
from numbers import Integral, Real
from pathlib import Path
from typing import Literal

import matplotlib as mpl
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.colors import LinearSegmentedColormap, is_color_like, to_hex
from plotnine import (
    element_blank,
    element_line,
    element_rect,
    element_text,
    theme,
    theme_classic,
    theme_set,
)
from plotnine.options import get_option, set_option

from ._compat import SKIP_GGANN_FRAMES
from .theme import _resolve_font, sizes, theme_ggann

__all__ = [
    "PublicationStyle",
    "publication_style",
    "theme_publication",
    "style_context",
    "publication_palette",
    "save_publication",
]


_QUALITATIVE = (
    "#E69F00",  # orange
    "#56B4E9",  # sky blue
    "#009E73",  # bluish green
    "#F0E442",  # yellow
    "#0072B2",  # blue
    "#D55E00",  # vermillion
    "#CC79A7",  # reddish purple
    "#2B2B2B",  # near black
)

# Extra colours are deterministic, but the first eight are the reviewed core
# vocabulary. Beyond eight categories users should add redundant encodings or
# review the palette in the figure's final context.
_QUALITATIVE_EXTENDED = _QUALITATIVE + (
    "#88CCEE",
    "#CC6677",
    "#DDCC77",
    "#117733",
    "#332288",
    "#AA4499",
    "#44AA99",
    "#999933",
    "#882255",
    "#661100",
    "#6699CC",
    "#AA4466",
)

_GOLDEN = 0.6180339887498949
_EXT_LIGHTNESS = (0.45, 0.62, 0.52)
_EXT_SATURATION = (0.70, 0.75, 0.55)


def _extension_colour(index: int) -> str:
    """The ``index``-th colour of the fixed sequence used past the vocabulary.

    A pure function of ``index`` (golden-ratio hue steps, cycling lightness), so a
    palette of ``n + 1`` colours always begins with the palette of ``n``.
    """
    hue = (_GOLDEN * (index + 1)) % 1.0
    lightness = _EXT_LIGHTNESS[index % 3]
    saturation = _EXT_SATURATION[index % 3]
    return to_hex(colorsys.hls_to_rgb(hue, lightness, saturation), keep_alpha=False)


def _extend_palette(base: Sequence[str], n: int) -> tuple[str, ...]:
    """Prefix-stable extension of ``base`` to ``n`` distinct colours."""
    colours = [to_hex(colour, keep_alpha=False) for colour in base]
    index = 0
    while len(colours) < n:
        candidate = _extension_colour(index)
        index += 1
        if candidate not in colours:
            colours.append(candidate)
    return tuple(colours[:n])


_FONT_FALLBACK = ("Arial", "Helvetica", "DejaVu Sans", "sans-serif")
_MISSING_COLOR = "#B3B3B3"


@dataclass(frozen=True, slots=True)
class PublicationStyle:
    """Validated, immutable journal-oriented design settings.

    Sizes are expressed in final-output units: physical dimensions in
    millimetres and typography/lines in points. The presets are generic layout
    conveniences, not claims of compliance with a particular journal.

    Parameters
    ----------
    preset : str
        Human-readable preset name stored in figure manifests.
    width_mm, height_mm : float
        Default final canvas dimensions in millimetres.
    base_size : float, default=6.5
        Base text size in points.
    axis_text_size, axis_title_size : float
        Tick-label and axis-title sizes in points.
    legend_text_size, legend_title_size : float
        Legend text and title sizes in points.
    strip_text_size, title_size, tag_size : float
        Facet-strip, plot-title, and panel-tag sizes in points.
    line_width : float, default=0.5
        Axis, tick, and guide line width in points.
    tick_length : float, default=2.0
        Major-tick length in points.
    point_size : float, default=1.2
        Recommended default point size for direct backends.
    legend_key_size : float, default=8.0
        Legend-key size in points.
    dpi : {300, 600}, default=600
        Raster resolution used by the style.
    fonts : tuple of str
        Ordered sans-serif fallback families. The default is Arial, Helvetica,
        DejaVu Sans, then the generic sans-serif family.
    missing_color : str
        Neutral colour for missing observations and matrix cells.
    qualitative : tuple of str
        Ordered categorical colour vocabulary.
    sequential_cmap : str
        Matplotlib colormap name used for unsigned continuous values.
    diverging : tuple of three str
        Low, midpoint, and high colours for zero-centred values.
    tag_levels : {"A", "a", "1", "i"}
        Default automatic panel-tag sequence.

    Returns
    -------
    PublicationStyle
        Frozen, validated style instance.

    Raises
    ------
    ValueError
        If a size is not positive, a colour or colormap is invalid, the DPI is
        unsupported, or the tag sequence is unknown.

    Notes
    -----
    Instances are frozen dataclasses. Use :func:`publication_style` to derive a
    modified copy rather than changing fields in place.

    Examples
    --------
    >>> style = PublicationStyle("custom", 120, 90, base_size=7)
    >>> style.width_mm
    120
    """

    preset: str
    width_mm: float
    height_mm: float
    base_size: float = 6.5
    axis_text_size: float = 5.5
    axis_title_size: float = 6.5
    legend_text_size: float = 5.5
    legend_title_size: float = 6.0
    strip_text_size: float = 6.2
    title_size: float = 7.0
    tag_size: float = 8.0
    line_width: float = 0.5
    tick_length: float = 2.0
    point_size: float = 1.2
    legend_key_size: float = 8.0
    dpi: int = 600
    fonts: tuple[str, ...] = _FONT_FALLBACK
    missing_color: str = _MISSING_COLOR
    qualitative: tuple[str, ...] = _QUALITATIVE_EXTENDED
    sequential_cmap: str = "viridis"
    diverging: tuple[str, str, str] = ("#2166AC", "#F7F7F7", "#B35806")
    tag_levels: str = "a"

    def __post_init__(self) -> None:
        if not isinstance(self.preset, str) or not self.preset:
            raise ValueError("preset must be a non-empty string.")
        for name in ("fonts", "qualitative", "diverging"):
            value = getattr(self, name)
            if isinstance(value, str) or not isinstance(value, Sequence):
                raise ValueError(f"{name} must be a sequence of strings.")
            object.__setattr__(self, name, tuple(value))
        positive = {
            "width_mm": self.width_mm,
            "height_mm": self.height_mm,
            "base_size": self.base_size,
            "axis_text_size": self.axis_text_size,
            "axis_title_size": self.axis_title_size,
            "legend_text_size": self.legend_text_size,
            "legend_title_size": self.legend_title_size,
            "strip_text_size": self.strip_text_size,
            "title_size": self.title_size,
            "tag_size": self.tag_size,
            "line_width": self.line_width,
            "tick_length": self.tick_length,
            "point_size": self.point_size,
            "legend_key_size": self.legend_key_size,
        }
        for name, value in positive.items():
            if (
                isinstance(value, bool)
                or not isinstance(value, Real)
                or not math.isfinite(float(value))
                or float(value) <= 0
            ):
                raise ValueError(f"{name} must be a positive finite number, got {value!r}.")
        if (
            isinstance(self.dpi, bool)
            or not isinstance(self.dpi, Integral)
            or self.dpi
            not in {
                300,
                600,
            }
        ):
            raise ValueError(f"dpi must be 300 or 600, got {self.dpi!r}.")
        if not self.fonts or not all(isinstance(font, str) and font for font in self.fonts):
            raise ValueError("fonts must contain at least one non-empty family name.")
        if not self.qualitative:
            raise ValueError("qualitative must contain at least one colour.")
        colours = [self.missing_color, *self.qualitative, *self.diverging]
        if not all(isinstance(colour, str) for colour in colours):
            raise ValueError("publication colours must be strings.")
        invalid = [colour for colour in colours if not is_color_like(colour)]
        if invalid:
            raise ValueError(f"Invalid publication colour(s): {invalid!r}.")
        if len(self.diverging) != 3:
            raise ValueError("diverging must contain low, midpoint, and high colours.")
        if self.tag_levels not in {"A", "a", "1", "i"}:
            raise ValueError("tag_levels must be 'A', 'a', '1', or 'i'.")
        if not isinstance(self.sequential_cmap, str):
            raise ValueError("sequential_cmap must be a Matplotlib colormap name.")
        try:
            mpl.colormaps[self.sequential_cmap]
        except KeyError as exc:
            raise ValueError(f"Unknown sequential_cmap {self.sequential_cmap!r}.") from exc

    def to_dict(self) -> dict[str, object]:
        """Return a serialisable copy suitable for manifests."""
        return asdict(self)


_PRESETS: Mapping[str, Mapping[str, object]] = {
    "single-column": {"preset": "single-column", "width_mm": 89.0, "height_mm": 70.0},
    "double-column": {"preset": "double-column", "width_mm": 183.0, "height_mm": 120.0},
}


def publication_style(
    preset: str | PublicationStyle = "single-column", **overrides: object
) -> PublicationStyle:
    """Create a validated publication style from a generic width preset.

    Parameters
    ----------
    preset : {"single-column", "double-column"} or PublicationStyle
        Starting style. Passing an existing style is useful for small immutable
        variations.
    **overrides
        Any :class:`PublicationStyle` field.

    Returns
    -------
    PublicationStyle
        New immutable style instance.

    Raises
    ------
    TypeError
        If an override name is unknown.
    ValueError
        If the preset or a field value is invalid.

    Notes
    -----
    ``"single-column"`` is 89 mm wide and ``"double-column"`` is 183 mm
    wide. These are ergonomic starting points; check the destination journal's
    current author instructions before submission.

    Examples
    --------
    >>> style = publication_style("double-column", base_size=7)
    >>> (style.width_mm, style.base_size)
    (183.0, 7)
    """
    if isinstance(preset, PublicationStyle):
        try:
            return replace(preset, **overrides)
        except TypeError as exc:
            raise TypeError(f"Unknown PublicationStyle override: {exc}") from exc
    if preset not in _PRESETS:
        choices = ", ".join(sorted(_PRESETS))
        raise ValueError(f"preset must be one of {choices}; got {preset!r}.")
    values = dict(_PRESETS[preset])
    values.update(overrides)
    try:
        return PublicationStyle(**values)
    except TypeError as exc:
        raise TypeError(f"Unknown PublicationStyle override: {exc}") from exc


_ACTIVE_STYLE: ContextVar[PublicationStyle | None] = ContextVar(
    "ggann_publication_style", default=None
)


def _active_style() -> PublicationStyle | None:
    return _ACTIVE_STYLE.get()


def _publication_active() -> bool:
    return _active_style() is not None


def _resolve_style(style: PublicationStyle | str | None) -> PublicationStyle:
    if style is None:
        return _active_style() or publication_style()
    return publication_style(style) if isinstance(style, str) else style


def theme_publication(
    style: PublicationStyle | str | None = None,
    *,
    axes: Literal["standard", "embedding", "matrix", "distribution"] = "standard",
):
    """Return a normal plotnine theme for final-size journal figures.

    Parameters
    ----------
    style : PublicationStyle or {"single-column", "double-column"}, optional
        Style instance or preset. The active context is used when omitted.
    axes : {"standard", "embedding", "matrix", "distribution"}
        Family-specific axis treatment. Embeddings suppress arbitrary units,
        matrices remove redundant lines, and distributions rotate x labels.

    Returns
    -------
    plotnine.theme
        A regular additive plotnine theme.

    Raises
    ------
    ValueError
        If ``axes`` or a style preset is invalid.

    Notes
    -----
    The theme is composable: a later user theme, scale, coordinate, annotation,
    layer, or facet remains authoritative.

    Examples
    --------
    >>> from plotnine import theme
    >>> final = exploratory + theme_publication(axes="embedding")
    >>> final = final + theme(legend_position="bottom")
    """
    if axes not in {"standard", "embedding", "matrix", "distribution"}:
        raise ValueError("axes must be 'standard', 'embedding', 'matrix', or 'distribution'.")
    resolved = _resolve_style(style)
    result = theme_classic(base_size=resolved.base_size, base_family="sans-serif")
    result += theme(
        text=element_text(family="sans-serif", color="#1A1A1A"),
        plot_background=element_rect(fill="white", color=None),
        panel_background=element_rect(fill="white", color=None),
        panel_grid=element_blank(),
        panel_border=element_blank(),
        axis_line=element_line(color="#1A1A1A", size=resolved.line_width),
        axis_ticks=element_line(color="#1A1A1A", size=resolved.line_width),
        axis_ticks_length=resolved.tick_length,
        axis_text=element_text(size=resolved.axis_text_size, color="#1A1A1A"),
        axis_title=element_text(size=resolved.axis_title_size, color="#1A1A1A"),
        legend_background=element_rect(fill="white", color=None),
        legend_key=element_blank(),
        legend_key_size=resolved.legend_key_size,
        legend_text=element_text(size=resolved.legend_text_size),
        legend_title=element_text(size=resolved.legend_title_size),
        strip_background=element_blank(),
        strip_text=element_text(size=resolved.strip_text_size, weight="bold"),
        plot_title=element_text(size=resolved.title_size, weight="bold", ha="left"),
        plot_tag=element_text(size=resolved.tag_size, weight="bold", ha="left"),
        plot_margin=0.018,
        figure_size=(resolved.width_mm / 25.4, resolved.height_mm / 25.4),
        dpi=resolved.dpi,
    )
    if axes == "embedding":
        result += theme(
            axis_title=element_blank(),
            axis_text=element_blank(),
            axis_ticks=element_blank(),
            axis_line=element_blank(),
        )
    elif axes == "matrix":
        result += theme(
            axis_title=element_blank(),
            axis_ticks=element_blank(),
            axis_line=element_blank(),
            axis_text_x=element_text(size=resolved.axis_text_size, rotation=45, ha="right"),
        )
    elif axes == "distribution":
        result += theme(
            axis_text_x=element_text(size=resolved.axis_text_size, rotation=30, ha="right")
        )
    # Keep this as the final margin override: multiline colourbar titles sit
    # above their bar and need more headroom in compact composition cells.
    result += theme(plot_margin_top=0.06)
    result._ggann_publication_style = resolved
    return result


def _family_theme(
    family: Literal["standard", "embedding", "matrix", "distribution"] = "standard",
):
    """Select the frozen legacy theme or the active publication family theme."""
    style = _active_style()
    return theme_ggann() if style is None else theme_publication(style, axes=family)


def _style_rcparams(style: PublicationStyle) -> dict[str, object]:
    return {
        "font.family": "sans-serif",
        "font.sans-serif": list(style.fonts),
        "font.size": style.base_size,
        "axes.labelsize": style.axis_title_size,
        "axes.titlesize": style.title_size,
        "axes.linewidth": style.line_width,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "xtick.labelsize": style.axis_text_size,
        "ytick.labelsize": style.axis_text_size,
        "xtick.major.width": style.line_width,
        "ytick.major.width": style.line_width,
        "xtick.major.size": style.tick_length,
        "ytick.major.size": style.tick_length,
        "legend.fontsize": style.legend_text_size,
        "legend.title_fontsize": style.legend_title_size,
        "legend.frameon": False,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "savefig.bbox": None,
    }


@contextmanager
def style_context(
    style: PublicationStyle | str | None = None, **overrides: object
) -> Iterator[PublicationStyle]:
    """Temporarily apply publication defaults and restore all global state.

    Parameters
    ----------
    style : PublicationStyle or {"single-column", "double-column"}, optional
        Base style. An omitted value inherits the enclosing context, if any.
    **overrides
        Field overrides used to derive the active immutable style.

    Returns
    -------
    contextlib.AbstractContextManager
        Context manager yielding the active :class:`PublicationStyle`.

    Raises
    ------
    TypeError
        If an override name is unknown.
    ValueError
        If a preset or override value is invalid.

    Notes
    -----
    Contexts may be nested. Plotnine's current theme, :data:`ggann.sizes`, every
    Matplotlib rcParam, and the active ggann style are restored after normal or
    exceptional exit. Plots retain the theme and direct-backend specification
    captured while they were built.

    Restoring is a full rollback to the state at entry: a :func:`set_theme` call
    or an ``rcParams`` edit made *inside* the block is discarded on exit, so make
    such changes before entering the block if they must persist.

    The active style is held in a ``contextvars.ContextVar``, which new threads do
    not inherit, whereas the plotnine theme, ``sizes`` and ``rcParams`` it swaps
    are process-global and visible to every thread while the block runs. Behaviour
    with concurrent threads is therefore unverified; use the context from one
    thread at a time.

    A requested font family that is not installed (for example Arial) triggers one
    ``UserWarning`` per process and text falls back to the next installed family.

    Examples
    --------
    >>> with style_context("double-column"):
    ...     panel = plot_embedding(adata, color="cell_type")
    """
    base = _active_style() if style is None else style
    resolved = publication_style(base or "single-column", **overrides)
    # Keep the raw option, including ``None`` or a theme class. Calling
    # ``theme_get`` would materialise a fresh default and would therefore not
    # restore the exact pre-context state.
    previous_theme = get_option("current_theme")
    previous_sizes = dict(vars(sizes))
    previous_rc = deepcopy(dict(mpl.rcParams))
    _resolve_font(resolved.fonts, context="style_context")
    token = _ACTIVE_STYLE.set(resolved)
    try:
        sizes.update(resolved.base_size)
        for key, value in _style_rcparams(resolved).items():
            mpl.rcParams[key] = value
        theme_set(theme_publication(resolved))
        yield resolved
    finally:
        set_option("current_theme", previous_theme)
        vars(sizes).clear()
        vars(sizes).update(previous_sizes)
        mpl.rcParams.update(previous_rc)
        _ACTIVE_STYLE.reset(token)


def publication_palette(
    kind: Literal["qualitative", "sequential", "diverging"],
    n: int | None = None,
    *,
    categories: Sequence[object] | None = None,
    style: PublicationStyle | str | None = None,
) -> tuple[str, ...] | dict[object, str]:
    """Return a deterministic publication colour vocabulary.

    Parameters
    ----------
    kind : {"qualitative", "sequential", "diverging"}
        Colour-vocabulary type.
    n : int, optional
        Number of colours. Defaults to eight for qualitative and 256 otherwise.
    categories : sequence, optional
        Ordered category values. Qualitative palettes then return a mapping
        instead of a tuple.
    style : PublicationStyle or {"single-column", "double-column"}, optional
        Style supplying the colour vocabulary.

    Returns
    -------
    tuple of str or dict
        Hex colours, or a stable category-to-colour mapping when ``categories``
        is supplied.

    Raises
    ------
    ValueError
        If ``kind`` is unknown, ``n`` is invalid, or categories are requested
        for a non-qualitative palette.

    Notes
    -----
    Missing observations use :attr:`PublicationStyle.missing_color` in plotting
    scales and are not added as a synthetic category. The core eight qualitative
    colours are accessibility-tested; larger vocabularies need visual review or
    a redundant encoding. Colours are assigned in the order of ``categories``, so
    reordering categories reassigns colours. Past the style's vocabulary the
    palette is extended from a fixed sequence, so ``n + 1`` colours always start
    with the ``n`` colours.

    Examples
    --------
    >>> colours = publication_palette(
    ...     "qualitative", categories=["B cell", "T cell", "NK cell"]
    ... )
    >>> colours["T cell"]
    '#56b4e9'
    """
    if kind not in {"qualitative", "sequential", "diverging"}:
        raise ValueError("kind must be 'qualitative', 'sequential', or 'diverging'.")
    resolved = _resolve_style(style)
    if categories is not None:
        if kind != "qualitative":
            raise ValueError("categories can only be supplied for a qualitative palette.")
        categories = list(categories)
        try:
            if len(set(categories)) != len(categories):
                raise ValueError("categories must contain unique values.")
        except TypeError as exc:
            raise ValueError("categories must contain hashable values.") from exc
        if n is not None and n != len(categories):
            raise ValueError("n must equal len(categories) when both are supplied.")
        n = len(categories)
    if n is None:
        n = 8 if kind == "qualitative" else 256
    if isinstance(n, bool) or not isinstance(n, Integral) or n < 1:
        raise ValueError(f"n must be a positive integer, got {n!r}.")
    n = int(n)

    if kind == "qualitative":
        # Past the vocabulary, extend from a fixed sequence so growing a category
        # list never recolours earlier categories. Separation is measured only for
        # the reviewed vocabulary; larger sets need a redundant encoding.
        colours = _extend_palette(resolved.qualitative[:n], n)
    elif kind == "sequential":
        cmap = mpl.colormaps[resolved.sequential_cmap].resampled(n)
        colours = tuple(to_hex(cmap(i), keep_alpha=False) for i in range(n))
    else:
        cmap = LinearSegmentedColormap.from_list("ggann_diverging", resolved.diverging, N=n)
        colours = tuple(to_hex(cmap(i), keep_alpha=False) for i in range(n))
    if categories is not None:
        return dict(zip(categories, colours, strict=True))
    return colours


_FORMATS = {"svg", "pdf", "png", "tiff"}
_RASTER_FORMATS = {"png", "tiff"}
_UNIT_TO_INCH = {"mm": 1 / 25.4, "cm": 1 / 2.54, "in": 1.0}


def _normalise_format(value: str) -> str:
    value = value.lower().lstrip(".")
    if value == "tif":
        value = "tiff"
    if value not in _FORMATS:
        raise ValueError("formats must contain only SVG, PDF, PNG, or TIFF.")
    return value


_KNOWN_SUFFIXES = _FORMATS | {"tif"}
_REJECTED_SUFFIXES = {"jpg", "jpeg", "eps", "ps", "gif", "bmp", "webp"}


def _output_paths(filename: str | Path, formats: Sequence[str] | None) -> list[tuple[str, Path]]:
    """Resolve output paths.

    Only a recognised format suffix (svg, pdf, png, tiff, tif) is treated as a
    suffix. Any other dot stays in the stem, so ``fig1.v2`` writes ``fig1.v2.svg``,
    ``fig1.v2.pdf`` and ``fig1.v2.png``. A format suffix that disagrees with an
    explicit ``formats`` raises instead of silently writing another format.
    """
    path = Path(filename)
    suffix = path.suffix.lower().lstrip(".")
    if suffix in _REJECTED_SUFFIXES:
        raise ValueError("formats must contain only SVG, PDF, PNG, or TIFF.")
    has_format_suffix = suffix in _KNOWN_SUFFIXES
    if formats is None:
        requested = [_normalise_format(suffix)] if has_format_suffix else ["svg", "pdf", "png"]
    else:
        requested = [_normalise_format(value) for value in formats]
        if not requested:
            raise ValueError("formats cannot be empty.")
        if has_format_suffix and _normalise_format(suffix) not in requested:
            raise ValueError(
                f"filename suffix '.{suffix}' disagrees with formats={tuple(formats)!r}; "
                "drop the suffix or include it in formats."
            )
    if len(set(requested)) != len(requested):
        raise ValueError("formats cannot contain duplicates.")
    stem = path.stem if has_format_suffix else path.name
    priority = {"svg": 0, "pdf": 1, "png": 2, "tiff": 3}
    requested.sort(key=priority.__getitem__)
    return [(fmt, path.with_name(f"{stem}.{fmt}")) for fmt in requested]


def _validate_dpi(dpi) -> int:
    """Export DPI: any finite positive integer (style presets stay 300/600)."""
    if isinstance(dpi, bool) or not isinstance(dpi, Integral) or dpi <= 0:
        raise ValueError(f"dpi must be a positive integer, got {dpi!r}.")
    return int(dpi)


def _width_in_units(
    width: str | float,
    units: Literal["mm", "cm", "in"],
) -> tuple[float, PublicationStyle]:
    if isinstance(width, str):
        style = publication_style(width)
        return style.width_mm / 25.4, style
    if (
        isinstance(width, bool)
        or not isinstance(width, Real)
        or not math.isfinite(float(width))
        or float(width) <= 0
    ):
        raise ValueError(f"width must be a positive number or preset, got {width!r}.")
    return float(width) * _UNIT_TO_INCH[units], _active_style() or publication_style()


def _captured_style(plot) -> PublicationStyle | None:
    """Find a publication style retained by a supported plotting object."""
    style = getattr(plot, "_ggann_publication_style", None)
    if isinstance(style, PublicationStyle):
        return style
    spec = getattr(plot, "_ggann_render_spec", None)
    style = getattr(spec, "publication_style", None)
    if isinstance(style, PublicationStyle):
        return style
    plot_theme = getattr(plot, "theme", None)
    style = getattr(plot_theme, "_ggann_publication_style", None)
    if isinstance(style, PublicationStyle):
        return style
    from plotnine.composition import Compose

    if isinstance(plot, Compose):
        for child in plot:
            if child_style := _captured_style(child):
                return child_style
    return None


def _figure_from_object(plot, width_in: float, height_in: float, dpi: int):
    """Draw a copied supported object at the requested canvas size."""
    from matplotlib.figure import Figure
    from plotnine import ggplot
    from plotnine import theme as p9_theme
    from plotnine.composition import Compose

    if isinstance(plot, Figure):
        figure = deepcopy(plot)
        figure.set_size_inches(width_in, height_in, forward=True)
        figure.set_dpi(dpi)
        return figure

    working = deepcopy(plot)
    if isinstance(working, Compose):
        working = working + p9_theme(figure_size=(width_in, height_in), dpi=dpi)
        # plotnine 0.15 creates the composition figure before entering any
        # child plot theme. Seed the canvas size here so physical gaps are
        # calculated against the requested export dimensions.
        with mpl.rc_context({"figure.figsize": (width_in, height_in), "figure.dpi": dpi}):
            figure = working.draw(show=False)
    elif isinstance(working, ggplot):
        # Direct Matplotlib ggplot subclasses deliberately invalidate their fast
        # path when grammar components are added. Let them draw their captured
        # specification, then resize the normalised axes canvas.
        if getattr(working, "fast_path_active", False):
            figure = working.draw(show=False)
        else:
            working += p9_theme(figure_size=(width_in, height_in), dpi=dpi)
            figure = working.draw(show=False)
    else:
        figure = None
        for name in ("figure", "fig"):
            candidate = getattr(working, name, None)
            if isinstance(candidate, Figure):
                figure = candidate
                break
        if figure is None:
            axis = getattr(working, "ax", None)
            candidate = getattr(axis, "figure", None)
            if isinstance(candidate, Figure):
                figure = candidate
        if figure is None and hasattr(working, "render"):
            working.render()
            for name in ("figure", "fig"):
                candidate = getattr(working, name, None)
                if isinstance(candidate, Figure):
                    figure = candidate
                    break
        if figure is None and callable(getattr(working, "plot", None)):
            import matplotlib.pyplot as plt

            figure = plt.figure(figsize=(width_in, height_in), dpi=dpi)
            working.plot()
        if figure is None:
            raise TypeError(
                "save_publication supports plotnine plots/compositions, Matplotlib "
                "figures, and ggann's documented grid-backend results."
            )
    figure.set_size_inches(width_in, height_in, forward=True)
    figure.set_dpi(dpi)
    if hasattr(working, "_apply_publication_gap"):
        figure.canvas.draw()
        working._apply_publication_gap(figure)
    for axis in list(getattr(figure, "_ggann_fit_left", ())):
        fit_left_labels(figure, axis)
    return figure


def fit_left_labels(figure, axis, *, pad_mm: float = 1.0, min_axes_fraction: float = 0.2) -> None:
    """Widen ``axis``'s left margin until its y tick labels fit inside the canvas.

    Labels are measured with the figure's renderer at its current physical size,
    not estimated from character counts. The right edge of the axes stays put and
    the margin only grows. The axis is remembered on the figure so
    :func:`save_publication` can repeat the fit after resizing the canvas.
    """
    renderer = figure._get_renderer()
    bbox = axis.yaxis.get_tightbbox(renderer)
    registered = figure.__dict__.setdefault("_ggann_fit_left", [])
    if axis not in registered:
        registered.append(axis)
    if bbox is None:
        return
    pad_px = pad_mm / 25.4 * figure.dpi
    need = (pad_px - bbox.x0) / figure.bbox.width
    if need <= 0:
        return
    pos = axis.get_position()
    width = max(pos.width - need, min_axes_fraction)
    axis.set_position((pos.x1 - width, pos.y0, width, pos.height))


def _overflow_report(figure, tolerance_px: float = 0.5) -> list[str]:
    """Describe text, legends and colourbars that leave the canvas or are clipped.

    Draws the figure once at its final physical size. Reports items beyond the
    canvas edge and texts cut by their own axes clip box. Nothing is cropped or
    resized.
    """
    from matplotlib.legend import Legend
    from matplotlib.offsetbox import AnchoredOffsetbox
    from matplotlib.text import Text

    # A fresh Agg canvas triggers plotnine's layout engine; redrawing the canvas a
    # plotnine figure already carries leaves artists at stale positions.
    FigureCanvasAgg(figure).draw()
    renderer = figure._get_renderer()
    canvas = figure.bbox
    mm_per_px = 25.4 / figure.dpi
    found: list[str] = []

    def beyond(box) -> str | None:
        sides = {
            "left": canvas.x0 - box.x0,
            "right": box.x1 - canvas.x1,
            "bottom": canvas.y0 - box.y0,
            "top": box.y1 - canvas.y1,
        }
        worst = max(sides, key=lambda side: sides[side])
        if sides[worst] <= tolerance_px:
            return None
        return f"{worst} by {sides[worst] * mm_per_px:.1f} mm"

    for artist in figure.findobj(lambda a: isinstance(a, (Text, Legend, AnchoredOffsetbox))):
        if not artist.get_visible():
            continue
        if isinstance(artist, Text):
            if not artist.get_text().strip():
                continue
            name = repr(artist.get_text().replace("\n", " ")[:30])
        else:
            name = "legend" if isinstance(artist, Legend) else "legend or colourbar"
        box = artist.get_window_extent(renderer)
        if not (box.width > 0 and box.height > 0):
            continue
        message = beyond(box)
        if message:
            found.append(f"{name} ({message})")
            continue
        if isinstance(artist, Text) and artist.get_clip_on():
            clip = artist.get_clip_box()
            if clip is not None and (
                box.x0 < clip.x0 - tolerance_px
                or box.x1 > clip.x1 + tolerance_px
                or box.y0 < clip.y0 - tolerance_px
                or box.y1 > clip.y1 + tolerance_px
            ):
                found.append(f"{name} (clipped by its axes)")
    return found


def _warn_overflow(figure, width_mm: float, height_mm: float) -> None:
    try:
        found = _overflow_report(figure)
    except Exception as exc:  # detection must never block an export
        warnings.warn(
            f"save_publication: overflow check skipped ({exc!r}).",
            skip_file_prefixes=SKIP_GGANN_FRAMES,
        )
        return
    if not found:
        return
    shown = "; ".join(found[:4]) + (f"; and {len(found) - 4} more" if len(found) > 4 else "")
    warnings.warn(
        f"save_publication: {len(found)} item(s) are cut off at the {width_mm:.0f} x "
        f"{height_mm:.0f} mm canvas: {shown}. The canvas was not changed; shorten labels, "
        "reduce font sizes or enlarge the figure.",
        UserWarning,
        skip_file_prefixes=SKIP_GGANN_FRAMES,
    )


def save_publication(
    plot,
    filename: str | Path,
    *,
    width: str | float = "single-column",
    height: float,
    units: Literal["mm", "cm", "in"] = "mm",
    formats: Sequence[str] | None = None,
    dpi: int = 600,
    background: str = "white",
) -> tuple[Path, ...]:
    """Save exact-size editable/vector and high-resolution raster outputs.

    Parameters
    ----------
    plot
        Plotnine plot/composition, Matplotlib figure, or a result returned by
        :func:`plot_clustermap` or :func:`plot_upset`.
    filename : str or pathlib.Path
        Output path or stem. Only an SVG, PDF, PNG, TIFF or TIF suffix counts as
        a format; any other dot stays in the stem, so ``fig1.v2`` writes
        ``fig1.v2.svg``, ``fig1.v2.pdf`` and ``fig1.v2.png``.
    width : float or {"single-column", "double-column"}
        Canvas width in ``units`` or a generic physical-width preset.
    height : float
        Positive canvas height in ``units``. Required: no default height is
        guessed from the width.
    units : {"mm", "cm", "in"}, default="mm"
        Units for numeric dimensions.
    formats : sequence of str, optional
        Any of SVG, PDF, PNG, or TIFF. A stem without a format suffix defaults to
        SVG, PDF, and PNG; a filename with a format suffix defaults to that one
        format. A suffix that is not in ``formats`` raises ``ValueError``.
    dpi : int, default=600
        Any finite positive integer resolution for PNG and TIFF; vector
        dimensions are physical. Raster size is ``floor(size_in_inches * dpi)``
        pixels (Matplotlib's rule), for example 89 mm at 150 dpi is 525 px wide.
        :class:`PublicationStyle` presets still accept only 300 or 600.
    background : str, default="white"
        ``"transparent"`` or any Matplotlib-compatible colour.

    Returns
    -------
    tuple of pathlib.Path
        Written paths in SVG, PDF, PNG, TIFF priority order.

    Raises
    ------
    TypeError
        If ``plot`` is not a supported plotting result.
    ValueError
        If dimensions, units, formats, DPI, or background are invalid, or a
        filename suffix disagrees with ``formats``.

    Warns
    -----
    UserWarning
        If text, a legend or a colourbar extends beyond the canvas or is clipped
        by its axes at the final size, or if the requested font family is not
        installed. The canvas is never enlarged or cropped.

    Notes
    -----
    The input is copied before rendering, the canvas is never tight-cropped,
    temporary figures are closed, SVG text remains text (``svg.fonttype="none"``),
    and PDF text uses embedded TrueType fonts (``pdf.fonttype=42``). A font family
    that is not installed falls back to the next installed one (a warning is
    issued once); the SVG still names the requested family, while the PDF embeds
    the font actually used. SVG and PDF files carry different metadata (the SVG
    stores a creation date and element ids, the PDF a producer string and
    timestamps), so byte-identical output across runs is not expected. A
    DPI-only change affects PNG and TIFF; SVG and PDF are resolution-independent.
    Axis-level guide collection from :func:`compose` stays in one panel.

    Examples
    --------
    >>> paths = save_publication(
    ...     figure, "figure_1", width="double-column", height=120, dpi=600
    ... )
    >>> [path.suffix for path in paths]
    ['.svg', '.pdf', '.png']
    """
    if units not in _UNIT_TO_INCH:
        raise ValueError("units must be 'mm', 'cm', or 'in'.")
    if (
        isinstance(height, bool)
        or not isinstance(height, Real)
        or not math.isfinite(float(height))
        or float(height) <= 0
    ):
        raise ValueError(f"height must be a positive number, got {height!r}.")
    dpi = _validate_dpi(dpi)
    if background != "transparent" and not is_color_like(background):
        raise ValueError(f"background must be white, transparent, or a colour; got {background!r}.")

    if isinstance(formats, str):
        raise ValueError("formats must be a sequence such as ('svg', 'pdf', 'png').")
    width_in, width_style = _width_in_units(width, units)
    style = _captured_style(plot) or _active_style() or width_style
    height_in = float(height) * _UNIT_TO_INCH[units]
    outputs = _output_paths(filename, formats)
    _resolve_font(style.fonts, context="save_publication")
    transparent = background == "transparent"
    facecolor = "none" if transparent else background
    figure = None
    try:
        with mpl.rc_context(_style_rcparams(style)):
            figure = _figure_from_object(plot, width_in, height_in, dpi)
            _warn_overflow(figure, width_in * 25.4, height_in * 25.4)
            saved: list[Path] = []
            for fmt, path in outputs:
                path.parent.mkdir(parents=True, exist_ok=True)
                figure.savefig(
                    path,
                    format=fmt,
                    dpi=dpi if fmt in _RASTER_FORMATS else None,
                    facecolor=facecolor,
                    edgecolor=facecolor,
                    transparent=transparent,
                    bbox_inches=None,
                    pad_inches=0,
                )
                saved.append(path)
            return tuple(saved)
    finally:
        if figure is not None:
            import matplotlib.pyplot as plt

            plt.close(figure)
