"""Tests for set_theme / reset_theme / sizes (the publication init)."""

from __future__ import annotations

import plotnine as p9
import pytest

import ggann as ag


@pytest.fixture(autouse=True)
def _restore_default_theme():
    # set_theme mutates plotnine's GLOBAL default; restore it after each test so
    # the change never leaks into other test files.
    yield
    ag.reset_theme()


def test_sizes_scale_and_geom_converter():
    ag.set_theme(base_size=10, register=False)
    assert ag.sizes.normal == 10.0
    assert ag.sizes.small == 8.0
    assert ag.sizes.large == 12.5
    assert ag.sizes.title == 13.0
    # pt -> mm for geom_text (ggplot2 .pt = 72.27/25.4)
    assert ag.sizes.geom(10) == pytest.approx(10 / 2.845276, rel=1e-6)


def test_set_theme_returns_theme_and_registers_default():
    ag.reset_theme()
    before = p9.theme_get()
    returned = ag.set_theme(base_size=9)
    assert isinstance(returned, p9.theme)
    # registered as the global default -> theme_get() changed
    assert p9.theme_get() is not before
    assert ag.sizes.normal == 9.0


def test_registered_default_applies_to_bare_ggplot():
    import pandas as pd

    ag.set_theme()
    df = pd.DataFrame({"x": [1, 2, 3], "y": [1, 2, 3]})
    (p9.ggplot(df, p9.aes("x", "y")) + p9.geom_point())._build()  # no explicit theme


def test_register_false_does_not_change_default():
    ag.reset_theme()
    before = p9.theme_get()
    ag.set_theme(base_size=12, register=False)
    assert p9.theme_get() is before  # default untouched
    assert ag.sizes.normal == 12.0  # but sizes still updated


def test_helpers_still_build_after_set_theme(adata, markers, group_key):
    ag.set_theme()
    (ag.plot_dotplot(adata, markers, group_key))._build()
    (ag.plot_violin(adata, markers[:1], group_key))._build()


def test_set_theme_base_size_reaches_helper_plots(adata, markers, group_key):
    import ggann as ag

    def axis_title_size():
        figure = ag.plot_dotplot(adata, markers, group_key).draw()
        return figure.axes[0].xaxis.get_label().get_fontsize()

    default = axis_title_size()
    try:
        ag.set_theme(base_size=16, register=False)
        assert axis_title_size() == 16.0 != default
    finally:
        ag.reset_theme()
    assert axis_title_size() == default
