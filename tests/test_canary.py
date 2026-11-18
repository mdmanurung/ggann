"""Canaries for the plotnine and plotnine-extra coupling.

Run alone with ``pytest -m canary``; they also run in the normal test job, and
plotnine-extra's ``Nightly canary`` workflow runs them against the newest
plotnine pre-release.
They assert effects, not imports, because the plotnine 0.16.0a15 breaks were
silent: ``compose(widths=...)`` stopped changing panel sizes and nothing raised.

A failure here points at ``src/ggann/_plotnine_compat.py``.
"""

from __future__ import annotations

from io import BytesIO

import anndata as ad
import numpy as np
import pandas as pd
import plotnine as p9
import plotnine_extra as pe
import pytest

import ggann as ag
from ggann import _plotnine_compat

pytestmark = pytest.mark.canary

# Names ggann re-exports from plotnine-extra (see ggann/__init__.py).
REEXPORTED = (
    "Beside",
    "Stack",
    "Wrap",
    "geom_label_repel",
    "geom_signif",
    "geom_text_repel",
    "ggsave",
    "plot_annotation",
    "plot_layout",
    "stat_anova_test",
    "stat_central_tendency",
    "stat_compare_means",
    "stat_cor",
    "stat_kruskal_test",
    "stat_pvalue_manual",
    "stat_pwc",
    "stat_regline_equation",
)


def _panel() -> p9.ggplot:
    frame = pd.DataFrame({"x": [1, 2, 3], "y": [1, 2, 3]})
    return p9.ggplot(frame, p9.aes("x", "y")) + p9.geom_point()


def _render(figure) -> None:
    """plotnine 0.16 lays panels out when the figure renders, not in draw()."""
    figure.savefig(BytesIO(), format="png")


def test_plotnine_hooks_ggann_relies_on_still_exist():
    assert _plotnine_compat.missing_coupling() == []


@pytest.mark.parametrize("name", REEXPORTED)
def test_reexport_is_the_plotnine_extra_object(name):
    assert getattr(ag, name) is getattr(pe, name)


def test_compose_widths_change_the_rendered_panels():
    fig = ag.compose([_panel(), _panel()], ncol=2, widths=[1, 3])
    _render(fig.draw(show=False))
    left, right = (p.axs[0].get_position().width for p in fig.iter_plots_all())
    assert right / left == pytest.approx(3, rel=0.05)


def test_compose_heights_change_the_rendered_panels():
    fig = ag.compose([_panel(), _panel()], ncol=1, heights=[1, 1.5])
    _render(fig.draw(show=False))
    top, bottom = (p.axs[0].get_position().height for p in fig.iter_plots_all())
    assert bottom / top == pytest.approx(1.5, rel=0.05)


def test_compose_gap_is_the_requested_millimetres():
    fig = ag.compose([_panel(), _panel()], ncol=2, gap=3.0)
    rendered = fig.draw(show=False)
    _render(rendered)
    spec = _plotnine_compat.sub_gridspec(fig)
    left, right = spec[0].get_position(rendered), spec[1].get_position(rendered)
    gap_mm = (right.x0 - left.x1) * rendered.get_figwidth() * 25.4
    assert gap_mm == pytest.approx(3.0, abs=0.01)


def test_gganndata_plot_draws_through_plotnine():
    rng = np.random.default_rng(0)
    obs = pd.DataFrame(
        {"group": pd.Categorical(np.tile(["a", "b"], 6))},
        index=[f"cell_{i}" for i in range(12)],
    )
    adata = ad.AnnData(rng.poisson(1.0, (12, 3)).astype("float32"), obs=obs)
    adata.var_names = ["g1", "g2", "g3"]
    adata.obsm["X_pca"] = rng.normal(size=(12, 2))
    mapping = ag.aes(x=ag.obsm("pca", 0), y=ag.obsm("pca", 1), color="group")
    plot = ag.gganndata(adata, mapping) + p9.geom_point()
    assert len(plot.layer_data(0)) == adata.n_obs
    _render(plot.draw(show=False))


# Layer class names the shapes module adds; each shapes PR appends its own.
SHAPE_LAYER_NAMES = ["geom_curve", "stat_curve"]


@pytest.mark.parametrize("name", SHAPE_LAYER_NAMES)
def test_plotnine_has_no_clashing_shape_names(name):
    """Guards ggann/shapes: plotnine's registry is keyed by class name, so a
    same-named class upstream would replace ggann's silently."""
    assert not hasattr(p9, name)
    assert not hasattr(pe, name)


def test_geom_curve_draws_one_bent_line_per_row():
    """Guards ggann/shapes/curve.py: the stat's ``group`` at draw, ``xend`` and
    ``yend`` on a path layer, a stat class in ``DEFAULT_PARAMS``, and kwargs
    reaching stat params."""
    df = pd.DataFrame({"x": [0, 0], "y": [0, 0], "xend": [4, 0], "yend": [0, 2]})
    plot = p9.ggplot(df, p9.aes("x", "y", xend="xend", yend="yend")) + ag.geom_curve(
        curvature=0.5, n=5
    )
    ax = plot.draw().axes[0]
    assert len(ax.lines) == 2
    assert [len(line.get_xdata()) for line in ax.lines] == [5, 5]
    assert min(ax.lines[0].get_ydata()) == pytest.approx(-1.0)
    assert min(ax.lines[1].get_ydata()) == pytest.approx(0.0)
