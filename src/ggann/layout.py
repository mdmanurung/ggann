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

import functools
import math
from typing import Sequence

from plotnine import labs

__all__ = ["compose", "tag_panels"]


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
        start = ord(levels)
        return [chr(start + i) for i in range(n)]  # A..Z / a..z (up to 26)
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
    levels : {"A", "a", "1", "I", "i"}
        Tag sequence style.

    Returns
    -------
    list
        New tagged plot objects in input order.

    Raises
    ------
    ValueError
        If ``levels`` is unsupported.

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
    tag_levels: str | None = "A",
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
    tag_levels : {"A", "a", "1", "I", "i"}, optional
        Tag style, or ``None`` to omit tags.

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
    if tag_levels:
        panels = tag_panels(panels, tag_levels)
    n = len(panels)
    if ncol is None:
        ncol = math.ceil(n / nrow) if nrow else math.ceil(math.sqrt(n))
    rows = [panels[i : i + ncol] for i in range(0, n, ncol)]
    row_objs = [functools.reduce(lambda a, b: a | b, r) for r in rows]
    return functools.reduce(lambda a, b: a / b, row_objs)
