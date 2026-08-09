"""Contracts for publication styles, palettes, context state, and export."""

from __future__ import annotations

import struct
from dataclasses import FrozenInstanceError

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import plotnine as p9
import pytest
from plotnine.options import get_option

import ggann as ag


def _plot():
    frame = pd.DataFrame({"x": [1, 2, 3], "y": [1, 4, 2], "group": ["a", "b", "a"]})
    return (
        p9.ggplot(frame, p9.aes("x", "y", color="group"))
        + p9.geom_point()
        + p9.labs(x="Time (days)", y="Response", color="Group")
    )


def _png_dimensions(path) -> tuple[int, int]:
    payload = path.read_bytes()
    assert payload[:8] == b"\x89PNG\r\n\x1a\n"
    return struct.unpack(">II", payload[16:24])


def test_publication_style_is_immutable_and_validated():
    style = ag.publication_style(
        "double-column",
        base_size=7,
        fonts=["DejaVu Sans", "sans-serif"],
        diverging=["#2166ac", "#f7f7f7", "#b35806"],
    )
    assert style.width_mm == 183
    assert style.base_size == 7
    assert style.fonts == ("DejaVu Sans", "sans-serif")
    assert isinstance(style.diverging, tuple)
    with pytest.raises(FrozenInstanceError):
        style.base_size = 9  # type: ignore[misc]
    with pytest.raises(ValueError, match="positive"):
        ag.publication_style(width_mm=0)
    with pytest.raises(ValueError, match="300 or 600"):
        ag.publication_style(dpi=150)
    with pytest.raises(ValueError, match="preset"):
        ag.publication_style("journal-x")
    for invalid in (True, float("inf"), float("nan")):
        with pytest.raises(ValueError, match="positive finite"):
            ag.publication_style(base_size=invalid)


def test_theme_publication_is_composable_and_user_theme_wins():
    plot = (
        _plot()
        + ag.theme_publication(axes="distribution")
        + p9.theme(axis_text_x=p9.element_text(size=13))
    )
    assert isinstance(plot.theme, p9.theme)
    assert plot.theme.themeables["axis_text_x"].theme_element.properties["size"] == 13


def test_style_context_restores_theme_sizes_rcparams_and_nested_state():
    raw_theme = get_option("current_theme")
    old_sizes = vars(ag.sizes).copy()
    old_rc = {key: mpl.rcParams[key] for key in ("svg.fonttype", "pdf.fonttype", "font.size")}

    with ag.style_context("double-column") as outer:
        assert outer.width_mm == 183
        assert mpl.rcParams["svg.fonttype"] == "none"
        assert mpl.rcParams["pdf.fonttype"] == 42
        with ag.style_context(base_size=5.8) as inner:
            assert inner.width_mm == 183
            assert ag.sizes.normal == 5.8
        assert ag.sizes.normal == outer.base_size

    assert get_option("current_theme") is raw_theme
    assert vars(ag.sizes) == old_sizes
    assert {key: mpl.rcParams[key] for key in old_rc} == old_rc


def test_style_context_restores_after_exception():
    raw_theme = get_option("current_theme")
    old_rc = mpl.rcParams["svg.fonttype"]
    with pytest.raises(RuntimeError, match="sentinel"):
        with ag.style_context():
            raise RuntimeError("sentinel")
    assert get_option("current_theme") is raw_theme
    assert mpl.rcParams["svg.fonttype"] == old_rc


def test_publication_palette_is_stable_and_category_aware():
    first = ag.publication_palette("qualitative", 8)
    second = ag.publication_palette("qualitative", categories=["B", "T", "NK"])
    assert len(first) == 8
    assert second == {"B": first[0], "T": first[1], "NK": first[2]}
    assert len(ag.publication_palette("sequential", 11)) == 11
    assert len(ag.publication_palette("diverging", 9)) == 9
    with pytest.raises(ValueError, match=r"len\(categories\)"):
        ag.publication_palette("qualitative", 2, categories=["a"])
    with pytest.raises(ValueError, match="unique"):
        ag.publication_palette("qualitative", categories=["a", "a"])
    with pytest.raises(ValueError, match="positive integer"):
        ag.publication_palette("qualitative", True)


def test_publication_obs_palette_validation_does_not_mutate_adata(adata, group_key):
    key = f"{group_key}_colors"
    original = list(adata.uns[key])
    adata.uns[key] = ["not-a-colour"]
    with ag.style_context():
        with pytest.warns(UserWarning, match="exactly one valid colour"):
            mapping = ag.obs_colors(adata, group_key)
    assert len(mapping) == len(adata.obs[group_key].cat.categories)
    assert adata.uns[key] == ["not-a-colour"]
    adata.uns[key] = original


def test_save_publication_exact_canvas_editable_text_and_no_figure_leak(tmp_path):
    before = set(plt.get_fignums())
    paths = ag.save_publication(
        _plot(),
        tmp_path / "figure",
        width="single-column",
        height=70,
        dpi=300,
    )
    assert [path.suffix for path in paths] == [".svg", ".pdf", ".png"]
    width, height = _png_dimensions(paths[-1])
    assert abs(width - round(89 / 25.4 * 300)) <= 1
    assert abs(height - round(70 / 25.4 * 300)) <= 1
    svg = paths[0].read_text(encoding="utf-8")
    assert "<text" in svg
    assert "<image" not in svg
    pdf = paths[1].read_bytes()
    assert b"/Type3" not in pdf
    assert b"/FontFile2" in pdf or b"/CIDFontType2" in pdf
    assert set(plt.get_fignums()) == before


def test_save_publication_suffix_and_validation(tmp_path):
    paths = ag.save_publication(_plot(), tmp_path / "single.tiff", height=70, width=89, dpi=600)
    assert paths == (tmp_path / "single.tiff",)
    with pytest.raises(ValueError, match="300 or 600"):
        ag.save_publication(_plot(), tmp_path / "bad.png", height=70, dpi=150)
    with pytest.raises(ValueError, match="background"):
        ag.save_publication(_plot(), tmp_path / "bad.png", height=70, background="not-a-colour")
    with pytest.raises(ValueError, match="sequence"):
        ag.save_publication(_plot(), tmp_path / "bad", height=70, formats="svg")
    with pytest.raises(ValueError, match="positive"):
        ag.save_publication(_plot(), tmp_path / "bad.png", height=True)


def test_export_uses_style_captured_by_theme_or_context(tmp_path, monkeypatch):
    import ggann.publication as publication

    style = ag.publication_style("double-column", base_size=7.25)
    explicit = _plot() + ag.theme_publication(style)
    with ag.style_context(style):
        contextual = _plot()

    observed: list[float] = []
    original = publication._figure_from_object

    def record(*args, **kwargs):
        observed.append(float(mpl.rcParams["font.size"]))
        return original(*args, **kwargs)

    monkeypatch.setattr(publication, "_figure_from_object", record)
    for index, plot in enumerate((explicit, contextual)):
        ag.save_publication(plot, tmp_path / f"captured-{index}.png", height=70, dpi=300)
    assert observed == [style.base_size, style.base_size]


def test_publication_helper_keeps_prepared_data_and_changes_only_style(adata, markers, group_key):
    legacy = ag.plot_dotplot(adata, markers, group_key)
    with ag.style_context():
        publication = ag.plot_dotplot(adata, markers, group_key)
    pd.testing.assert_frame_equal(legacy.data, publication.data, check_exact=True)
    assert legacy.theme.themeables["text"].theme_element.properties["size"] == 11
    assert publication.theme.themeables["text"].theme_element.properties["size"] == 6.5


def test_embedding_publication_defaults_and_missing_values_are_represented(adata, group_key):
    selected = adata[:30].copy()
    selected.obs.iloc[0, selected.obs.columns.get_loc(group_key)] = None
    with ag.style_context():
        plot = ag.plot_embedding(selected, color=group_key, pointdensity=False)
    assert len(plot.data) == selected.n_obs
    assert plot.data[group_key].isna().sum() == 1
    assert plot.layers[0].geom._kwargs["raster"] is False


def test_direct_backend_captures_publication_style_and_rasterization(tmp_path, adata, group_key):
    with ag.style_context():
        plot = ag.plot_embedding(
            adata,
            color=group_key,
            pointdensity=False,
            backend="matplotlib",
            rasterized=True,
        )
    assert plot.fast_path_active
    assert plot._ggann_render_spec.publication_style.base_size == 6.5
    path = ag.save_publication(
        plot,
        tmp_path / "rasterized.svg",
        width="single-column",
        height=70,
        dpi=300,
    )[0]
    svg = path.read_text(encoding="utf-8")
    assert svg.count("<image") == 1
    assert "<text" in svg


def test_tile_rasterization_retains_vector_annotations(tmp_path, adata, markers, group_key):
    with ag.style_context():
        plot = ag.plot_matrixplot(
            adata,
            markers,
            group_key,
            rasterized=True,
            annotate="force",
        )
    path = ag.save_publication(
        plot,
        tmp_path / "matrix.svg",
        width="single-column",
        height=70,
        dpi=300,
    )[0]
    svg = path.read_text(encoding="utf-8")
    assert svg.count("<image") == 1
    assert "<text" in svg


def test_auto_annotations_require_twelve_point_cells(tmp_path):
    from ggann._annotation import geom_contrast_text

    frame = pd.DataFrame(
        [
            {"x": x, "y": y, "value": (x + y) / 50, "label": "CELLANNOTATION"}
            for x in range(25)
            for y in range(25)
        ]
    )
    base = (
        p9.ggplot(frame, p9.aes("x", "y", fill="value"))
        + p9.geom_tile()
        + p9.scale_fill_cmap(cmap_name="viridis")
    )
    automatic = base + geom_contrast_text(p9.aes(label="label"), min_cell_pt=12, format_string="{}")
    forced = base + geom_contrast_text(p9.aes(label="label"), min_cell_pt=0, format_string="{}")
    auto_path = ag.save_publication(automatic, tmp_path / "auto.svg", width=40, height=40, dpi=300)[
        0
    ]
    force_path = ag.save_publication(forced, tmp_path / "force.svg", width=40, height=40, dpi=300)[
        0
    ]
    assert "CELLANNOTATION" not in auto_path.read_text(encoding="utf-8")
    assert "CELLANNOTATION" in force_path.read_text(encoding="utf-8")


def test_grid_backend_exceptions_keep_return_types_and_export(tmp_path, adata, markers, group_key):
    from marsilea.upset import Upset
    from PyComplexHeatmap import ClusterMapPlotter

    with ag.style_context():
        clustermap = ag.plot_clustermap(
            adata,
            markers[:3],
            group_by=group_key,
            plot=False,
        )
        upset = ag.plot_upset(
            {"T cell": {"CD3D", "IL7R"}, "NK cell": {"NKG7", "GNLY"}},
            render=False,
        )
    assert isinstance(clustermap, ClusterMapPlotter)
    assert isinstance(upset, Upset)
    assert clustermap._ggann_publication_style.base_size == 6.5
    assert upset._ggann_publication_style.base_size == 6.5
    cm_path = ag.save_publication(
        clustermap, tmp_path / "clustermap.svg", width=89, height=70, dpi=300
    )[0]
    upset_path = ag.save_publication(upset, tmp_path / "upset.svg", width=89, height=70, dpi=300)[0]
    assert cm_path.stat().st_size > 0
    assert upset_path.stat().st_size > 0
