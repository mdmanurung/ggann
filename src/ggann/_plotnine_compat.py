"""The one place ggann reaches into private plotnine names.

Everything underscore-prefixed that ggann needs from plotnine lives here, so a
plotnine bump breaks this file and ``tests/test_canary.py`` names what moved.

Prefer plotnine's public API. Add to this module only when no public route
exists, and say why. ``missing_coupling()`` lists any hook or attribute below
that plotnine no longer has. Several of these fail silently when they go
missing (an overridden hook is simply never called), so the canary tests also
check the *effect*, for example that ``compose(widths=...)`` changes panel
sizes.
"""

from __future__ import annotations

from dataclasses import replace
from typing import Any

from plotnine import coord_radial, ggplot
from plotnine._mpl.gridspec import p9GridSpec
from plotnine.composition import Compose
from plotnine.mapping._env import Environment
from plotnine.stats.stat import stat

__all__ = [
    "BUILD_CACHE_ATTRS",
    "Environment",
    "figure_edge_fractions",
    "grid_shape",
    "missing_coupling",
    "pin_grid_spacing",
    "stat",
    "sub_gridspec",
]

# ``stat`` is the stat base class. plotnine exports neither ``plotnine.stat`` nor
# the class from ``plotnine.stats`` (that name is the module), so the shape
# layers in ``ggann.shapes`` import it from here.

# Attributes plotnine stores on a plot after building it. ``MatplotlibGGPlot``
# shares them by reference when it copies a plot. ``_build_objs`` is a
# deprecated alias of ``built`` that plotnine still sets.
BUILD_CACHE_ATTRS: tuple[str, ...] = ("built", "_build_objs")


def sub_gridspec(composition: Compose) -> p9GridSpec:
    """Return the grid that holds a composition's items.

    plotnine 0.16 builds it in ``Compose._generate_gridspecs`` and keeps it in
    ``_sub_gridspec``. Before 0.16.0a15 it was ``gridspec``.
    """
    return composition._sub_gridspec  # type: ignore[attr-defined,no-any-return]


def grid_shape(composition: Compose) -> tuple[int, int]:
    """Return ``(nrow, ncol)`` of a composition's grid."""
    layout = composition.layout
    return int(layout.nrow), int(layout.ncol)  # type: ignore[arg-type]


def figure_edge_fractions(figure: Any) -> tuple[float, float]:
    """Return the share of figure width and height left inside the margins.

    plotnine reserves the root theme's ``plot_margin_*`` after a composition
    sets up its grid. Margins are fractions of the figure width; top and bottom
    scale by width / height.
    """
    width_left = height_left = 1.0
    theme = getattr(getattr(figure, "_owner", None), "theme", None)
    if theme is not None:
        fig_w, fig_h = theme.getp("figure_size")
        width_left -= theme.getp("plot_margin_left") + theme.getp("plot_margin_right")
        height_left -= (theme.getp("plot_margin_top") + theme.getp("plot_margin_bottom")) * (
            fig_w / fig_h
        )
    return width_left, height_left


def pin_grid_spacing(gridspec: p9GridSpec, **spacing: float) -> None:
    """Set ``wspace`` / ``hspace`` on a grid and keep them through rendering.

    ``p9GridSpec`` lost ``update()`` in 0.16.0a15, and the layout engine resets
    its spacing to 0 when the figure renders. So set the stored params now and
    re-apply the values to the params the engine passes in.
    """
    gridspec._subplot_params.update(**spacing)
    apply_params = gridspec.update_params_and_artists

    def keep_spacing(gsparams: Any) -> None:
        apply_params(replace(gsparams, **spacing))

    gridspec.update_params_and_artists = keep_spacing  # type: ignore[method-assign]


def missing_coupling() -> list[str]:
    """Return the plotnine attributes ggann relies on that are gone."""
    wanted = {
        "Compose._generate_gridspecs": hasattr(Compose, "_generate_gridspecs"),
        "Compose._sub_gridspec (annotation)": "_sub_gridspec"
        in getattr(Compose, "__annotations__", {}),
        "p9GridSpec.update_params_and_artists": hasattr(p9GridSpec, "update_params_and_artists"),
        "p9GridSpec.bbox_relative": hasattr(p9GridSpec, "bbox_relative"),
        "ggplot.layer_data": hasattr(ggplot, "layer_data"),
        "ggplot.build": hasattr(ggplot, "build"),
        "stat.compute_panel": hasattr(stat, "compute_panel"),
        "coord_radial.munch": hasattr(coord_radial, "munch"),
        "coord_radial.is_linear": hasattr(coord_radial, "is_linear"),
        "Environment.capture": hasattr(Environment, "capture"),
    }
    return [name for name, present in wanted.items() if not present]
