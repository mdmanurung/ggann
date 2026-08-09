"""Deterministic final-size rendering and resource-lifetime contracts."""

from __future__ import annotations

import gc
import itertools
import tracemalloc
import weakref

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import plotnine as p9
import pytest
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.collections import PathCollection, PolyCollection
from matplotlib.figure import Figure
from matplotlib.text import Text
from matplotlib.transforms import Bbox
from PIL import Image

import ggann as ag
from ggann.publication import _figure_from_object, _style_rcparams


def _contract_plot(style):
    frame = pd.DataFrame(
        {
            "x": [1, 2, 3],
            "y": [1, 4, 2],
            "lineage": ["T cell", "B cell", "T cell"],
        }
    )
    return (
        p9.ggplot(frame, p9.aes("x", "y", color="lineage"))
        + p9.geom_point()
        + p9.labs(
            title="Measured response",
            tag=style.tag_levels,
            x="Time (days)",
            y="Response",
            color="Lineage",
        )
        + ag.theme_publication(style)
    )


def _visible_text(figure: Figure) -> list[tuple[Text, Bbox]]:
    canvas = FigureCanvasAgg(figure)
    canvas.draw()
    renderer = canvas.get_renderer()
    return [
        (artist, artist.get_window_extent(renderer))
        for artist in figure.findobj(Text)
        if artist.get_visible() and artist.get_text().strip()
    ]


@pytest.mark.parametrize(
    ("preset", "width_mm", "height_mm"),
    [("single-column", 89, 70), ("double-column", 183, 120)],
)
def test_final_size_text_bounds_minimum_size_and_overlap(preset, width_mm, height_mm):
    style = ag.publication_style(preset)
    plot = _contract_plot(style)
    with mpl.rc_context(_style_rcparams(style)):
        figure = _figure_from_object(plot, width_mm / 25.4, height_mm / 25.4, 144)
    try:
        texts = _visible_text(figure)
        tolerance = 0.5 * 144 / 72  # half a point in display pixels
        for artist, bounds in texts:
            assert artist.get_fontsize() >= 5
            assert bounds.x0 >= -tolerance
            assert bounds.y0 >= -tolerance
            assert bounds.x1 <= figure.bbox.x1 + tolerance
            assert bounds.y1 <= figure.bbox.y1 + tolerance

        for (_, first), (_, second) in itertools.combinations(texts, 2):
            width = max(0.0, min(first.x1, second.x1) - max(first.x0, second.x0))
            height = max(0.0, min(first.y1, second.y1) - max(first.y0, second.y0))
            smaller = min(first.width * first.height, second.width * second.height)
            if smaller:
                assert width * height / smaller < 0.05
    finally:
        plt.close(figure)


@pytest.mark.parametrize(
    ("preset", "width_mm", "height_mm"),
    [("single-column", 89, 70), ("double-column", 183, 120)],
)
def test_publication_pngs_are_exact_and_deterministic(tmp_path, preset, width_mm, height_mm):
    style = ag.publication_style(preset)
    plot = _contract_plot(style)
    first = ag.save_publication(
        plot,
        tmp_path / f"{preset}-first.png",
        width=preset,
        height=height_mm,
        dpi=300,
    )[0]
    second = ag.save_publication(
        plot,
        tmp_path / f"{preset}-second.png",
        width=preset,
        height=height_mm,
        dpi=300,
    )[0]
    with Image.open(first) as image:
        expected = (
            round(width_mm / 25.4 * 300),
            round(height_mm / 25.4 * 300),
        )
        assert all(abs(actual - target) <= 1 for actual, target in zip(image.size, expected))
    assert first.read_bytes() == second.read_bytes()


def test_tiff_has_exact_physical_dimensions(tmp_path):
    output = ag.save_publication(
        _contract_plot(ag.publication_style()),
        tmp_path / "figure.tiff",
        width="single-column",
        height=70,
        dpi=600,
    )[0]
    with Image.open(output) as image:
        expected = (round(89 / 25.4 * 600), round(70 / 25.4 * 600))
        assert all(abs(actual - target) <= 1 for actual, target in zip(image.size, expected))


def test_weighted_gap_composition_keeps_all_text_inside_canvas():
    style = ag.publication_style("double-column")
    composition = ag.compose(
        [_contract_plot(style) for _ in range(4)],
        ncol=2,
        widths=(0.9, 1.1),
        heights=(1, 1),
        gap=2,
        tag_levels="a",
    )
    with mpl.rc_context(_style_rcparams(style)):
        figure = _figure_from_object(composition, 183 / 25.4, 120 / 25.4, 144)
    try:
        texts = _visible_text(figure)
        tolerance = 0.5 * 144 / 72
        for _, bounds in texts:
            assert bounds.x0 >= -tolerance
            assert bounds.y0 >= -tolerance
            assert bounds.x1 <= figure.bbox.x1 + tolerance
            assert bounds.y1 <= figure.bbox.y1 + tolerance
    finally:
        plt.close(figure)


def test_rendered_point_and_tile_populations_match_prepared_data(adata, markers, group_key):
    selected = adata[:40].copy()
    with ag.style_context():
        embedding = ag.plot_embedding(
            selected,
            color=group_key,
            pointdensity=False,
        )
        matrix = ag.plot_matrixplot(selected, markers[:3], group_key)

    embedding_figure = embedding.draw(show=False)
    matrix_figure = matrix.draw(show=False)
    try:
        point_counts = [
            len(artist.get_offsets()) for artist in embedding_figure.findobj(PathCollection)
        ]
        tile_counts = [len(artist.get_paths()) for artist in matrix_figure.findobj(PolyCollection)]
        assert max(point_counts) == selected.n_obs == len(embedding.data)
        assert len(matrix.data) in tile_counts
    finally:
        plt.close(embedding_figure)
        plt.close(matrix_figure)


def test_twenty_exports_do_not_leak_figures_or_one_mib(tmp_path, monkeypatch):
    import ggann.publication as publication

    plot = _contract_plot(ag.publication_style())
    ag.save_publication(plot, tmp_path / "warm-1.png", height=70, dpi=300)
    ag.save_publication(plot, tmp_path / "warm-2.png", height=70, dpi=300)
    before_figures = set(plt.get_fignums())
    references: list[weakref.ReferenceType[Figure]] = []
    original = publication._figure_from_object

    def record(*args, **kwargs):
        figure = original(*args, **kwargs)
        references.append(weakref.ref(figure))
        return figure

    monkeypatch.setattr(publication, "_figure_from_object", record)
    gc.collect()
    tracemalloc.start()
    start = tracemalloc.get_traced_memory()[0]
    try:
        for _ in range(20):
            ag.save_publication(plot, tmp_path / "repeat.png", height=70, dpi=300)
        gc.collect()
        growth = tracemalloc.get_traced_memory()[0] - start
    finally:
        tracemalloc.stop()

    # A subsequent render releases Matplotlib/plotnine's last-render cache, so
    # every figure produced by the measured exports must then be collectible.
    monkeypatch.setattr(publication, "_figure_from_object", original)
    ag.save_publication(plot, tmp_path / "release-cache.png", height=70, dpi=300)
    gc.collect()
    assert all(reference() is None for reference in references)
    assert growth < 1024 * 1024
    assert set(plt.get_fignums()) == before_figures


def test_exception_path_closes_and_releases_temporary_figure(tmp_path, monkeypatch):
    import ggann.publication as publication

    references: list[weakref.ReferenceType[Figure]] = []
    original_build = publication._figure_from_object
    original_save = Figure.savefig

    def record(*args, **kwargs):
        figure = original_build(*args, **kwargs)
        references.append(weakref.ref(figure))
        return figure

    def fail_save(*args, **kwargs):
        raise RuntimeError("sentinel export failure")

    monkeypatch.setattr(publication, "_figure_from_object", record)
    monkeypatch.setattr(Figure, "savefig", fail_save)
    before = set(plt.get_fignums())
    with pytest.raises(RuntimeError, match="sentinel"):
        ag.save_publication(_contract_plot(ag.publication_style()), tmp_path / "bad.png", height=70)
    assert set(plt.get_fignums()) == before

    monkeypatch.setattr(Figure, "savefig", original_save)
    ag.save_publication(
        _contract_plot(ag.publication_style()),
        tmp_path / "next.png",
        height=70,
    )
    gc.collect()
    assert references[0]() is None
