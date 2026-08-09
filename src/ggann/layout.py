"""Assemble multi-panel figures with panel tags (A, B, C ...).

Borrows exactplot's figure-assembly idea. Uses plotnine's native composition
operators (``|`` side-by-side, ``/`` stacked) plus ``labs(tag=)`` for panel
labels, so it stays vector-clean with no extra dependency. Save at an exact
physical size with plotnine's own ``.save(width=, height=, units="mm")`` -- that
covers exactplot's millimetre workflow without tikz/LaTeX.

    fig = ag.compose([p_umap, p_dotplot, p_violin, p_props], ncol=2)
    fig.save("figure1.pdf", width=180, height=140, units="mm")
"""

from __future__ import annotations

import math
from copy import deepcopy
from dataclasses import dataclass
from io import BytesIO
from pathlib import Path
from typing import Sequence, cast

from plotnine import labs, theme
from plotnine.composition import Beside, Compose, Stack

__all__ = ["compose", "tag_panels"]


def _validate_ratios(name: str, values: Sequence[float] | None, length: int) -> tuple[float, ...]:
    if values is None:
        return (1.0,) * length
    values = tuple(values)
    if len(values) != length:
        raise ValueError(f"{name} must contain exactly {length} values, got {len(values)}.")
    if any(
        isinstance(value, bool)
        or not isinstance(value, (int, float))
        or not math.isfinite(float(value))
        or float(value) <= 0
        for value in values
    ):
        raise ValueError(f"{name} must contain only positive finite numbers.")
    return tuple(float(value) for value in values)


class _GGAnnComposition:
    """Shared exact-save and millimetre-gap behavior for native compositions."""

    _ggann_gap_mm: float

    def _configure_physical_gap(self, figure) -> None:
        """Set the parent grid gap before child plot layout is calculated."""
        if self._ggann_gap_mm <= 0:
            return
        gap_in = self._ggann_gap_mm / 25.4
        bbox = self.gridspec.bbox_relative
        updates: dict[str, float] = {}
        if self.ncol > 1:
            available = bbox.width * figure.get_figwidth()
            average_cell = (available - gap_in * (self.ncol - 1)) / self.ncol
            if average_cell <= 0:
                raise ValueError("gap is too large for the requested figure width.")
            updates["wspace"] = gap_in / average_cell
        if self.nrow > 1:
            available = bbox.height * figure.get_figheight()
            average_cell = (available - gap_in * (self.nrow - 1)) / self.nrow
            if average_cell <= 0:
                raise ValueError("gap is too large for the requested figure height.")
            updates["hspace"] = gap_in / average_cell
        self.gridspec.update(**updates)

    def _apply_publication_gap(self, figure) -> None:
        """Compatibility hook; gaps are already encoded in parent gridspecs."""
        figure._ggann_gap_applied = True

    def draw(self, *, show: bool = False):
        figure = super().draw(show=show)
        self._apply_publication_gap(figure)
        return figure

    def save(
        self,
        filename: str | Path | BytesIO,
        format: str | None = None,
        dpi: int | None = None,
        *,
        width: float | None = None,
        height: float | None = None,
        units: str = "in",
        **kwargs,
    ) -> None:
        """Save, honoring explicit dimensions unlike native plotnine 0.15."""
        if width is None and height is None:
            return super().save(filename, format=format, dpi=dpi, **kwargs)
        if width is None or height is None:
            raise ValueError("width and height must be supplied together.")
        if isinstance(filename, BytesIO):
            factors = {"in": 1.0, "cm": 1 / 2.54, "mm": 1 / 25.4}
            if units not in factors:
                raise ValueError("units must be 'in', 'cm', or 'mm'.")
            plot = cast(Compose, self) + theme(
                figure_size=(width * factors[units], height * factors[units]),
                dpi=dpi or 300,
            )
            figure = plot.draw(show=False)
            try:
                figure.savefig(filename, format=format, dpi=dpi or 300, bbox_inches=None)
            finally:
                import matplotlib.pyplot as plt

                plt.close(figure)
            return None

        from .publication import save_publication

        requested_formats = [format] if format else None
        save_publication(
            self,
            filename,
            width=width,
            height=height,
            units=units,
            formats=requested_formats,
            dpi=dpi or 300,
            background=(
                "transparent" if kwargs.get("transparent") else kwargs.get("facecolor", "white")
            ),
        )
        return None


@dataclass(repr=False)
class _WeightedBeside(_GGAnnComposition, Beside):
    _ggann_ratios: tuple[float, ...]
    _ggann_gap_mm: float

    def _create_gridspec(self, figure, nest_into):
        from plotnine._mpl.gridspec import p9GridSpec

        self.gridspec = p9GridSpec(
            self.nrow,
            self.ncol,
            figure,
            width_ratios=self._ggann_ratios,
            nest_into=nest_into,
        )
        self._configure_physical_gap(figure)


@dataclass(repr=False)
class _WeightedStack(_GGAnnComposition, Stack):
    _ggann_ratios: tuple[float, ...]
    _ggann_gap_mm: float

    def _create_gridspec(self, figure, nest_into):
        from plotnine._mpl.gridspec import p9GridSpec

        self.gridspec = p9GridSpec(
            self.nrow,
            self.ncol,
            figure,
            height_ratios=self._ggann_ratios,
            nest_into=nest_into,
        )
        self._configure_physical_gap(figure)


def _guide_signature(plot) -> tuple | None:
    """Return a trained guide contract, or ``None`` when it cannot be proven."""
    aesthetics = {"color", "colour", "fill", "size", "shape", "linetype", "alpha"}
    mappings: list[tuple[str, str]] = []
    for mapping in [plot.mapping, *(layer.mapping for layer in plot.layers)]:
        mappings.extend(
            sorted((name, str(value)) for name, value in mapping.items() if name in aesthetics)
        )
    if not mappings:
        return None

    working = deepcopy(plot)
    try:
        working._build()
    except Exception:
        # Guide collection is conservative: a plot that cannot be trained here
        # keeps its guide rather than risking a semantically incorrect merge.
        return None

    scale_signatures: list[tuple] = []
    for scale in working._build_objs.scales:
        used = tuple(sorted(set(scale.aesthetics) & aesthetics))
        if not used or scale.guide is None or scale.guide is False or scale.guide == "none":
            continue
        try:
            breaks = scale.get_breaks()
            labels = scale.get_labels(breaks)
            palette = scale.map(breaks)
        except Exception:
            return None
        scale_signatures.append(
            (
                type(scale).__module__,
                type(scale).__qualname__,
                used,
                repr(scale.final_limits),
                tuple(map(repr, breaks)),
                tuple(map(repr, labels)),
                tuple(map(repr, palette)),
                repr(scale.name),
                repr(scale.na_value),
                repr(scale.guide),
            )
        )
    if not scale_signatures:
        return None
    labels = tuple(
        sorted(
            (name, str(value)) for name, value in vars(working.labels).items() if name in aesthetics
        )
    )
    return tuple(mappings), tuple(scale_signatures), labels


def _collect_duplicate_guides(panels: list) -> list:
    seen: set[tuple] = set()
    collected = []
    for panel in panels:
        if isinstance(panel, Compose):
            # Native composition guide internals are not declarative enough to
            # compare safely; leave nested compositions untouched.
            collected.append(panel)
            continue
        signature = _guide_signature(panel)
        if signature is not None and signature in seen:
            panel = panel + theme(legend_position="none")
        elif signature is not None:
            seen.add(signature)
        collected.append(panel)
    return collected


def _roman(n: int) -> str:
    numerals = [(10, "x"), (9, "ix"), (5, "v"), (4, "iv"), (1, "i")]
    out = ""
    for value, sym in numerals:
        while n >= value:
            out += sym
            n -= value
    return out


def _tag_labels(levels: str, n: int) -> list[str]:
    if levels in ("A", "a"):
        if n > 26:
            raise ValueError(
                "alphabetic panel tags support at most 26 panels; "
                "use tag_levels='1' or 'i' for larger compositions"
            )
        start = ord(levels)
        return [chr(start + i) for i in range(n)]
    if levels == "1":
        return [str(i + 1) for i in range(n)]
    if levels == "i":
        return [_roman(i + 1) for i in range(n)]
    raise ValueError(f"tag_levels must be 'A', 'a', '1' or 'i', got {levels!r}")


def tag_panels(panels: Sequence, levels: str = "A") -> list:
    """Add panel tags to a list of plots.

    Parameters
    ----------
    panels : sequence
        Plotnine plots or compositions.
    levels : {"A", "a", "1", "i"}
        Tag sequence style.

    Returns
    -------
    list
        New tagged plot objects in input order.

    Raises
    ------
    ValueError
        If ``levels`` is unsupported, or alphabetic tags are requested for
        more than 26 panels.

    Notes
    -----
    Input plot objects are composed immutably; no AnnData is involved.

    Examples
    --------
    >>> tagged = tag_panels([p1, p2], levels="A")
    """
    panels = list(panels)
    tags = _tag_labels(levels, len(panels))
    return [p + labs(tag=t) for p, t in zip(panels, tags)]


def compose(
    panels: Sequence,
    *,
    ncol: int | None = None,
    nrow: int | None = None,
    tag_levels: str | None = "auto",
    widths: Sequence[float] | None = None,
    heights: Sequence[float] | None = None,
    gap: float = 0.0,
    guides: str = "keep",
):
    """Arrange plots into a tagged multi-panel figure.

    ``panels`` is a flat list of plotnine plots (ggann helpers return these). They
    are wrapped into a grid -- ``ncol`` / ``nrow`` control the shape (default: a
    roughly square layout) -- and tagged ``A``, ``B``, ... unless ``tag_levels=None``.

    Returns a plotnine composition; save it at an exact size with
    ``.save(width=, height=, units="mm")``. For uneven panel sizes, compose the
    sub-figures yourself with ``|`` / ``/`` and pass them in.

    Parameters
    ----------
    panels : sequence
        Plotnine plots or compositions.
    ncol, nrow : int, optional
        Grid dimensions; at most one is normally needed.
    tag_levels : {"auto", "A", "a", "1", "i"}, optional
        Tag style, or ``None`` to omit tags. ``"auto"`` keeps uppercase legacy
        tags and uses bold lowercase tags in a publication context.
    widths, heights : sequence of float, optional
        Positive relative column and row sizes.
    gap : float, default=0
        Physical gap between adjacent panels in millimetres.
    guides : {"keep", "collect"}, default="keep"
        Keep all guides or suppress only exact declarative duplicates.

    Returns
    -------
    plotnine.composition
        Composable multi-panel figure.

    Raises
    ------
    ValueError
        If panels are empty, dimensions are invalid, or tags are unsupported.

    Notes
    -----
    Composition does not render plots or mutate their source data.

    Examples
    --------
    >>> figure = compose([p1, p2], ncol=2)
    """
    panels = list(panels)
    if not panels:
        raise ValueError("compose needs at least one panel.")
    if guides not in {"keep", "collect"}:
        raise ValueError("guides must be 'keep' or 'collect'.")
    if (
        isinstance(gap, bool)
        or not isinstance(gap, (int, float))
        or not math.isfinite(float(gap))
        or float(gap) < 0
    ):
        raise ValueError("gap must be a non-negative finite number of millimetres.")
    if tag_levels == "auto":
        from .publication import _active_style

        style = _active_style()
        tag_levels = style.tag_levels if style is not None else "A"
    if tag_levels:
        panels = tag_panels(panels, tag_levels)
    if guides == "collect":
        panels = _collect_duplicate_guides(panels)
    n = len(panels)
    if ncol is not None and (not isinstance(ncol, int) or isinstance(ncol, bool) or ncol < 1):
        raise ValueError("ncol must be a positive integer.")
    if nrow is not None and (not isinstance(nrow, int) or isinstance(nrow, bool) or nrow < 1):
        raise ValueError("nrow must be a positive integer.")
    if ncol is not None and nrow is not None and ncol * nrow < n:
        raise ValueError("ncol * nrow is smaller than the number of panels.")
    if ncol is None:
        ncol = math.ceil(n / nrow) if nrow else math.ceil(math.sqrt(n))
    actual_nrow = math.ceil(n / ncol)
    if nrow is not None and actual_nrow > nrow:
        raise ValueError("nrow is too small for the requested ncol and panels.")
    width_ratios = _validate_ratios("widths", widths, ncol)
    height_ratios = _validate_ratios("heights", heights, actual_nrow)
    rows = [panels[i : i + ncol] for i in range(0, n, ncol)]
    row_objs = [
        row[0] if len(row) == 1 else _WeightedBeside(row, width_ratios[: len(row)], float(gap))
        for row in rows
    ]
    if len(row_objs) == 1:
        return row_objs[0]
    return _WeightedStack(row_objs, height_ratios, float(gap))
