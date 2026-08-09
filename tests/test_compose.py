"""Tests for compose / tag_panels (figure assembly)."""

from __future__ import annotations

import pandas as pd
import plotnine as p9
import pytest
from plotnine.composition import Compose

import ggann as ag
from ggann.layout import _tag_labels


def _panels(adata, group_key, markers, n=4):
    builders = [
        lambda: ag.plot_embedding(adata, "umap", color=group_key),
        lambda: ag.plot_dotplot(adata, markers, group_key),
        lambda: ag.plot_violin(adata, markers[:1], group_key),
        lambda: ag.plot_proportions(adata, group_key, split_by="phase"),
    ]
    return [b() for b in builders[:n]]


def test_tag_labels():
    assert _tag_labels("A", 3) == ["A", "B", "C"]
    assert _tag_labels("a", 2) == ["a", "b"]
    assert _tag_labels("1", 3) == ["1", "2", "3"]
    assert _tag_labels("i", 4) == ["i", "ii", "iii", "iv"]
    with pytest.raises(ValueError, match="at most 26"):
        _tag_labels("a", 27)
    with pytest.raises(ValueError, match="tag_levels"):
        _tag_labels("X", 2)


def test_compose_grid_builds(adata, group_key, markers):
    fig = ag.compose(_panels(adata, group_key, markers, 4), ncol=2)
    fig.save  # composition object is saveable
    # rendering the whole composition exercises every panel + the layout
    fig.draw(show=False)


def test_compose_default_and_single(adata, group_key, markers):
    ag.compose(_panels(adata, group_key, markers, 3)).draw(show=False)  # default shape
    ag.compose(_panels(adata, group_key, markers, 1)).draw(show=False)  # single panel


def test_compose_no_tags(adata, group_key, markers):
    ag.compose(_panels(adata, group_key, markers, 2), tag_levels=None).draw(show=False)


def test_tag_panels_returns_tagged_plots(adata, group_key, markers):
    tagged = ag.tag_panels(_panels(adata, group_key, markers, 2), levels="a")
    assert len(tagged) == 2
    assert all(isinstance(p, p9.ggplot) for p in tagged)


def test_compose_empty_raises():
    with pytest.raises(ValueError, match="at least one"):
        ag.compose([])


def test_compose_weighted_layout_is_native_composition(adata, group_key, markers):
    fig = ag.compose(
        _panels(adata, group_key, markers, 4),
        ncol=2,
        widths=[2, 1],
        heights=[1, 1.5],
        gap=1.5,
    )
    assert isinstance(fig, Compose)
    assert fig._ggann_ratios == (1.0, 1.5)
    rendered = fig.draw(show=False)
    panel_heights = [spec.plot.axs[0].get_position().height for spec in fig.plotspecs]
    assert sum(panel_heights[2:]) > sum(panel_heights[:2])
    assert rendered._ggann_gap_applied is True
    first_row = fig[0]
    left = first_row.gridspec[0].get_position(rendered)
    right = first_row.gridspec[1].get_position(rendered)
    horizontal_gap_mm = (right.x0 - left.x1) * rendered.get_figwidth() * 25.4
    top = fig.gridspec[0].get_position(rendered)
    bottom = fig.gridspec[1].get_position(rendered)
    vertical_gap_mm = (top.y0 - bottom.y1) * rendered.get_figheight() * 25.4
    assert horizontal_gap_mm == pytest.approx(1.5)
    assert vertical_gap_mm == pytest.approx(1.5)


def test_compose_layout_validation(adata, group_key, markers):
    panels = _panels(adata, group_key, markers, 2)
    with pytest.raises(ValueError, match="widths"):
        ag.compose(panels, ncol=2, widths=[1])
    with pytest.raises(ValueError, match="positive"):
        ag.compose(panels, ncol=2, widths=[1, 0])
    with pytest.raises(ValueError, match="gap"):
        ag.compose(panels, gap=-1)
    with pytest.raises(ValueError, match="guides"):
        ag.compose(panels, guides="merge")


def test_compose_collects_only_identical_guides():
    frame = pd.DataFrame({"x": [1, 2], "y": [2, 1], "group": ["a", "b"]})
    base = p9.ggplot(frame, p9.aes("x", "y", color="group")) + p9.geom_point()
    same = base + p9.scale_color_manual(values={"a": "red", "b": "blue"})
    different = base + p9.scale_color_manual(values={"a": "black", "b": "grey"})

    collected = ag.compose([same, same], ncol=2, guides="collect", tag_levels=None)
    assert collected[1].theme.getp("legend_position") == "none"
    kept = ag.compose([same, different], ncol=2, guides="collect", tag_levels=None)
    assert kept[1].theme.getp("legend_position") != "none"


def test_compose_trains_implicit_guides_before_collecting():
    first_data = pd.DataFrame({"x": [1, 2], "y": [2, 1], "group": ["a", "b"]})
    second_data = pd.DataFrame({"x": [1, 2], "y": [2, 1], "group": ["a", "c"]})
    first = p9.ggplot(first_data, p9.aes("x", "y", color="group")) + p9.geom_point()
    duplicate = p9.ggplot(first_data, p9.aes("x", "y", color="group")) + p9.geom_point()
    different = p9.ggplot(second_data, p9.aes("x", "y", color="group")) + p9.geom_point()

    collected = ag.compose([first, duplicate], ncol=2, guides="collect", tag_levels=None)
    assert collected[1].theme.getp("legend_position") == "none"
    kept = ag.compose([first, different], ncol=2, guides="collect", tag_levels=None)
    assert kept[1].theme.getp("legend_position") != "none"


def test_compose_save_honours_requested_mm_dimensions(tmp_path, adata, group_key, markers):
    import struct

    fig = ag.compose(_panels(adata, group_key, markers, 2), ncol=2, tag_levels=None)
    path = tmp_path / "composition.png"
    fig.save(path, width=183, height=120, units="mm", dpi=300)
    payload = path.read_bytes()
    width, height = struct.unpack(">II", payload[16:24])
    assert abs(width - round(183 / 25.4 * 300)) <= 1
    assert abs(height - round(120 / 25.4 * 300)) <= 1


def test_publication_compose_defaults_to_lowercase_tags(adata, group_key, markers):
    with ag.style_context():
        fig = ag.compose(_panels(adata, group_key, markers, 2), ncol=2)
    assert fig[0].labels.tag == "a"
    assert fig[1].labels.tag == "b"
