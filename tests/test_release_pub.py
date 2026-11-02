"""Regressions for the 0.1.0 readiness `pub` rows (theme, publication, layout, palette).

Row IDs refer to plans/evaluation/remaining_todos.md.
"""

from __future__ import annotations

import itertools
import re
import warnings

import matplotlib as mpl
import numpy as np
import pandas as pd
import plotnine as p9
import pytest
from PIL import Image

import ggann as ag
from ggann import _palette, theme
from ggann._compat import _resolve_alias
from ggann.publication import _output_paths, _overflow_report, _validate_dpi


@pytest.fixture(autouse=True)
def _restore_theme():
    yield
    ag.reset_theme()


def _points(**labels):
    frame = pd.DataFrame({"x": [1, 2, 3], "y": [1, 4, 2]})
    return p9.ggplot(frame, p9.aes("x", "y")) + p9.geom_point() + p9.labs(**labels)


# --------------------------------------------------------------------- B9 / S1 / S2
@pytest.mark.parametrize("bad", [0, -3, float("nan"), float("inf"), "9", None, True])
def test_B9_set_theme_rejects_bad_sizes_without_touching_state(bad):
    ag.set_theme(base_size=9)
    registered = p9.theme_get()
    with pytest.raises(ValueError, match="base_size"):
        ag.set_theme(base_size=bad)
    with pytest.raises(ValueError, match="base_size"):
        ag.theme_ggann(base_size=bad)
    assert ag.sizes.normal == 9.0
    assert p9.theme_get() is registered


def test_S1_reset_theme_resets_sizes():
    ag.set_theme(base_size=9)
    ag.reset_theme()
    assert ag.sizes.normal == 11.0
    assert ag.sizes.small == 8.8


def test_B9_dotted_stem_and_format_suffix_policy(tmp_path):
    stem = [p.name for _, p in _output_paths(tmp_path / "fig1.v2", None)]
    assert stem == ["fig1.v2.svg", "fig1.v2.pdf", "fig1.v2.png"]
    assert [p.name for _, p in _output_paths(tmp_path / "fig1.v2.PNG", None)] == ["fig1.v2.png"]
    with pytest.raises(ValueError, match="SVG, PDF, PNG, or TIFF"):
        _output_paths(tmp_path / "fig.jpg", None)


def test_S7_suffix_must_agree_with_formats(tmp_path):
    with pytest.raises(ValueError, match="disagrees with formats"):
        ag.save_publication(_points(), tmp_path / "x.png", height=40, formats=("svg",))
    assert not list(tmp_path.iterdir())
    ok = ag.save_publication(_points(), tmp_path / "x.png", height=40, formats=("svg", "png"))
    assert [p.name for p in ok] == ["x.svg", "x.png"]
    with pytest.raises(TypeError, match="height"):
        ag.save_publication(_points(), tmp_path / "y.png")  # type: ignore[call-arg]


def test_S2_style_context_rolls_back_in_block_changes():
    before_theme = p9.theme_get()
    before_size = ag.sizes.normal
    before_rc = mpl.rcParams["axes.linewidth"]
    with ag.style_context("single-column"):
        ag.set_theme(base_size=14)
        mpl.rcParams["axes.linewidth"] = 3.0
    assert p9.theme_get() is before_theme
    assert ag.sizes.normal == before_size
    assert mpl.rcParams["axes.linewidth"] == before_rc


# ---------------------------------------------------------------------- B5
def test_B5_numeric_obs_gets_continuous_scale_and_draws(adata):
    from plotnine.scales.scale_continuous import scale_continuous
    from plotnine.scales.scale_discrete import scale_discrete

    scale = ag.scale_color_obs(adata, "n_counts")
    assert isinstance(scale, scale_continuous)
    frame = adata.obs[["n_counts"]].assign(x=np.arange(adata.n_obs), y=1.0)
    p9.ggplot(frame, p9.aes("x", "y", color="n_counts")).__add__(p9.geom_point()).__add__(
        scale
    ).draw()

    flag = adata.copy()
    flag.obs["flag"] = np.arange(flag.n_obs) % 2 == 0
    assert isinstance(ag.scale_fill_obs(flag, "flag"), scale_discrete)
    assert isinstance(ag.scale_fill_obs(adata, "bulk_labels"), scale_discrete)
    with pytest.raises(KeyError, match="nope"):
        ag.scale_color_obs(adata, "nope")
    assert ag.obs_colors(adata, "n_counts") is None


# ---------------------------------------------------------------- D4 / S3
def _top_colour(scale):
    return mpl.colors.to_hex(scale.palette(np.array([1.0]))[0])


def test_S3_explicit_cmap_wins_inside_and_outside_style_context():
    reds_top = mpl.colors.to_hex(mpl.colormaps["Reds"](1.0))
    viridis_top = mpl.colors.to_hex(mpl.colormaps["viridis"](1.0))
    assert reds_top != viridis_top
    with ag.style_context("single-column"):
        explicit = ag.scale_color_expression("Reds")
        sentinel = ag.scale_color_expression()
        fill = ag.scale_fill_expression("viridis")
        fill_default = ag.scale_fill_expression()
    assert _top_colour(explicit) == reds_top
    assert _top_colour(sentinel) == viridis_top
    assert _top_colour(fill) == viridis_top == _top_colour(fill_default)
    assert _top_colour(ag.scale_color_expression()) == reds_top
    assert _top_colour(ag.scale_fill_expression()) == viridis_top


def test_D4_dotplot_cmap_precedence_and_signed_diverging(adata):
    genes = ["CD3D", "NKG7"]

    def scale(**kw):
        plot = ag.plot_dotplot(adata, genes, group_by="bulk_labels", use_raw=True, **kw)
        plot.draw()
        return next(s for s in plot.scales if "color" in s.aesthetics)

    reds_top = mpl.colors.to_hex(mpl.colormaps["Reds"](1.0))
    with ag.style_context("single-column"):
        assert _top_colour(scale(cmap="Reds")) == reds_top
        assert _top_colour(scale()) != reds_top
        assert type(scale(standard_scale="zscore")).__name__ == "scale_color_gradient2"
    assert _top_colour(scale()) == reds_top


# ----------------------------------------------------------- D3 / D8 / Q8
def test_D3_fallback_is_scanpy_default_and_stored_colours_win(adata):
    import warnings as w

    from scanpy.plotting import palettes

    plain = adata.copy()
    plain.uns.pop("bulk_labels_colors", None)
    cats = list(plain.obs["bulk_labels"].cat.categories)
    scale = ag.scale_fill_obs(plain, "bulk_labels")
    assert type(scale).__name__ == "scale_fill_manual"
    mapping = _palette._obs_mapping(plain, plain.obs["bulk_labels"], "bulk_labels")
    assert mapping == dict(zip(cats, palettes.default_20[: len(cats)]))

    stored = adata.copy()
    stored.uns["bulk_labels_colors"] = ["#%02x0000" % (10 * i + 20) for i in range(len(cats))]
    uns_before = list(stored.uns["bulk_labels_colors"])
    with w.catch_warnings():
        w.simplefilter("error")
        mapping = ag.obs_colors(stored, "bulk_labels")
    assert mapping is not None and list(mapping.values()) == uns_before
    assert stored.uns["bulk_labels_colors"] == uns_before


def test_D3_clustermap_annotation_reuses_stored_colours(adata):
    stored = adata.copy()
    cats = list(stored.obs["bulk_labels"].cat.categories)
    stored.uns["bulk_labels_colors"] = ["#%02x0000" % (10 * i + 20) for i in range(len(cats))]
    result = ag.plot_clustermap(
        stored, ["CD3D", "NKG7", "CST3"], group_by="bulk_labels", use_raw=True, plot=False
    )
    colours = result.top_annotation.annotations[0].color_dict
    assert colours == dict(zip(cats, stored.uns["bulk_labels_colors"]))


def test_D8_palette_extension_never_recolours_earlier_categories():
    prefixes = [ag.publication_palette("qualitative", n) for n in (8, 20, 21, 22, 40)]
    for shorter, longer in itertools.pairwise(prefixes):
        assert longer[: len(shorter)] == shorter
    assert len(set(prefixes[-1])) == 40
    categories = [f"c{i}" for i in range(22)]
    mapping = ag.publication_palette("qualitative", categories=categories)
    assert list(mapping.values()) == list(ag.publication_palette("qualitative", 22))


def test_Q8_duplicate_invalid_and_weak_stored_palettes_are_reported(adata):
    cats = list(adata.obs["bulk_labels"].cat.categories)
    dup = adata.copy()
    dup.uns["bulk_labels_colors"] = ["#ff0000"] * len(cats)
    _palette._WARNED.clear()
    with pytest.warns(UserWarning, match="duplicate colours"):
        mapping = ag.obs_colors(dup, "bulk_labels")
    assert set(mapping.values()) == {"#ff0000"}

    bad = adata.copy()
    bad.uns["bulk_labels_colors"] = ["notacolour"] * len(cats)
    with pytest.warns(UserWarning, match="invalid colours"):
        assert ag.obs_colors(bad, "bulk_labels") is None

    _palette._WARNED.clear()
    with ag.style_context("single-column"):
        with pytest.warns(UserWarning, match="colour-vision separation"):
            ag.obs_colors(adata, "bulk_labels")  # 10 categories: beyond the measured 8
        small = adata.copy()
        small.obs["g8"] = pd.Categorical(
            np.arange(small.n_obs) % 8, categories=list(range(8))
        ).astype(str)
        small.obs["g8"] = small.obs["g8"].astype("category")
        _palette._WARNED.clear()
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            ag.obs_colors(small, "g8")


# ------------------------------------------------------------------- S5
@pytest.mark.parametrize("bad", [True, 0, -150, float("nan"), float("inf"), 150.5, "150"])
def test_S5a_export_dpi_must_be_a_positive_integer(bad):
    with pytest.raises(ValueError, match="positive integer"):
        _validate_dpi(bad)


def test_S5a_compose_save_accepts_150_dpi_with_exact_pixels(tmp_path):
    figure = ag.compose([_points(), _points()], ncol=2, tag_levels=None)
    out = tmp_path / "draft.png"
    figure.save(out, width=89, height=60, units="mm", dpi=150)
    with Image.open(out) as image:
        assert image.size == (int(89 / 25.4 * 150), int(60 / 25.4 * 150))
    with pytest.raises(ValueError, match="positive integer"):
        figure.save(out, width=89, height=60, units="mm", dpi=0)
    with pytest.raises(ValueError, match="300 or 600"):
        ag.publication_style(dpi=150)


def test_S5c_compose_and_tag_panels_reject_non_plots():
    with pytest.raises(TypeError, match="panel 1 is a DataFrame"):
        ag.compose([_points(), pd.DataFrame()])
    with pytest.raises(TypeError, match="panel 0 is a str"):
        ag.tag_panels(["not a plot"])


# ------------------------------------------------------------------ S4a / S4b
def _drawn(plot, width_mm, height_mm, dpi=150):
    from ggann.publication import _figure_from_object

    return _figure_from_object(plot, width_mm / 25.4, height_mm / 25.4, dpi)


def test_S4a_overflow_is_reported_not_fixed(tmp_path):
    import matplotlib.pyplot as plt

    long = "A very long axis title that cannot fit on a 30 mm figure"
    tight = _points(x=long)
    figure = _drawn(tight, 30, 30)
    try:
        found = _overflow_report(figure)
    finally:
        plt.close(figure)
    assert any("very long" in item for item in found)

    with pytest.warns(UserWarning, match="cut off at the 30 x 30 mm canvas"):
        (path,) = ag.save_publication(tight, tmp_path / "t.png", width=30, height=30, dpi=100)
    with Image.open(path) as image:
        assert image.size == (int(30 / 25.4 * 100), int(30 / 25.4 * 100))  # canvas unchanged

    figure = _drawn(_points(x="x"), 89, 70)
    try:
        assert _overflow_report(figure) == []
    finally:
        plt.close(figure)
    with warnings.catch_warnings():
        warnings.filterwarnings("error", message=".*cut off", category=UserWarning)
        ag.save_publication(_points(x="x"), tmp_path / "ok.png", width=89, height=70, dpi=100)


def test_S4b_dotplot_labels_are_measured_at_70_mm(adata, monkeypatch):
    import matplotlib.pyplot as plt

    import ggann._matplotlib_backend as backend

    labelled = adata.copy()
    labelled.obs["ct"] = pd.Categorical(
        labelled.obs["bulk_labels"]
        .astype(str)
        .replace({"CD14+ Monocyte": "GR3A+ Monocytes (long label)"})
    )

    def left_clipped():
        with ag.style_context("single-column"):
            plot = ag.plot_dotplot(
                labelled, ["CD3D", "NKG7"], group_by="ct", use_raw=True, backend="matplotlib"
            )
            figure = _drawn(plot, 70, 70)
        try:
            return [item for item in _overflow_report(figure) if "left" in item]
        finally:
            plt.close(figure)

    assert left_clipped() == []
    # Without the measured fit (the old character-count margin) labels are cut off.
    monkeypatch.setattr(backend, "fit_left_labels", lambda *a, **k: None)
    monkeypatch.setattr("ggann.publication.fit_left_labels", lambda *a, **k: None)
    assert any("Monocytes" in item for item in left_clipped())


# ------------------------------------------------------------------- S6 / Q11
def test_S6_font_fallback_warns_once_and_unknown_family_is_reported(monkeypatch):
    monkeypatch.setattr(theme, "_WARNED_FONTS", set())
    theme._font_available.cache_clear()
    with pytest.warns(UserWarning, match="'NoSuchFont' is not installed"):
        ag.set_theme(family="NoSuchFont", register=False)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        ag.set_theme(family="NoSuchFont", register=False)  # once per process
        ag.set_theme(family="DejaVu Sans", register=False)
    monkeypatch.setattr(theme, "_font_available", lambda family: family != "Arial")
    with pytest.warns(
        UserWarning, match="'Arial' is not installed; text falls back to 'Helvetica'"
    ):
        with ag.style_context("single-column"):
            pass


def test_S6_svg_text_stays_editable_and_pdf_fonts_embedded(tmp_path):
    plot = _points(title="Hello") + ag.theme_publication()
    svg, pdf = ag.save_publication(
        plot, tmp_path / "f", width=89, height=70, formats=("svg", "pdf")
    )
    text = svg.read_text()
    assert text.count("<text") >= 3 and "font-family: 'Arial'" in text
    data = pdf.read_bytes()
    assert b"/Type3" not in data and re.search(rb"/FontFile2|/CIDFontType2", data)


# ---------------------------------------------------------------------- A1
def test_VERIFY4_grid_backends_export_and_clustermap_cut_off_is_reported(tmp_path, adata):
    with ag.style_context():
        clustermap = ag.plot_clustermap(
            adata, ["CD3D", "NKG7", "CST3"], group_by="bulk_labels", use_raw=True, plot=False
        )
        upset = ag.plot_upset({"T": {"CD3D", "IL7R"}, "NK": {"NKG7", "GNLY"}}, render=False)
    # an unplotted clustermap is plotted at export size, so it fits without a warning
    with warnings.catch_warnings():
        warnings.filterwarnings("error", message=".*cut off", category=UserWarning)
        ag.save_publication(clustermap, tmp_path / "cm.svg", width=89, height=70)
    assert (tmp_path / "cm.svg").stat().st_size > 0
    ag.save_publication(upset, tmp_path / "up.svg", width=89, height=70)
    assert (tmp_path / "up.svg").stat().st_size > 0


def _short_group_adata():
    import anndata as ad

    names = [f"T{i}" for i in range(6)]
    rng = np.random.default_rng(0)
    adata = ad.AnnData(
        rng.random((60, 6)).astype("float32"),
        obs=pd.DataFrame({"g": pd.Categorical(np.repeat(names, 10), categories=names)}),
        var=pd.DataFrame(index=[f"G{i}" for i in range(6)]),
    )
    adata.obs_names = [f"c{i}" for i in range(60)]
    return adata


@pytest.mark.parametrize(("case", "width", "height"), [("short", 89, 70), ("pbmc", 183, 120)])
def test_clustermap_autofit_keeps_everything_inside_the_exact_canvas(
    tmp_path, adata, case, width, height
):
    import matplotlib.pyplot as plt

    plt.figure(figsize=(width / 25.4, height / 25.4))  # PyComplexHeatmap lays out at this size
    if case == "short":
        data, genes, group = _short_group_adata(), [f"G{i}" for i in range(6)], "g"
    else:
        data, group = adata, "bulk_labels"
        genes = ["CD3D", "NKG7", "CST3", "GNLY", "MS4A1", "FCGR3A"]
    clustermap = ag.plot_clustermap(
        data,
        genes,
        group_by=group,
        use_raw=case == "pbmc",
        row_dendrogram=True,
        col_dendrogram=True,
    )
    try:
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "error", message=".*cut off", category=UserWarning
            )  # S4a: nothing is cut off
            (path,) = ag.save_publication(
                clustermap, tmp_path / "cm.png", width=width, height=height, dpi=150
            )
    finally:
        plt.close("all")
    with Image.open(path) as image:
        assert image.size == (int(width / 25.4 * 150), int(height / 25.4 * 150))


def test_clustermap_autofit_warns_when_the_labels_cannot_fit(tmp_path):
    import anndata as ad
    import matplotlib.pyplot as plt

    names = [f"a very long group name number {i}" for i in range(40)]
    rng = np.random.default_rng(0)
    adata = ad.AnnData(
        rng.random((80, 5)).astype("float32"),
        obs=pd.DataFrame({"g": pd.Categorical(np.repeat(names, 2), categories=names)}),
        var=pd.DataFrame(index=[f"G{i}" for i in range(5)]),
    )
    adata.obs_names = [f"c{i}" for i in range(80)]
    plt.figure(figsize=(30 / 25.4, 30 / 25.4))
    clustermap = ag.plot_clustermap(adata, list(adata.var_names), group_by="g")
    try:
        with pytest.warns(UserWarning, match="cut off at the 30 x 30 mm.*minimum body size"):
            (path,) = ag.save_publication(
                clustermap, tmp_path / "x.png", width=30, height=30, dpi=100
            )
    finally:
        plt.close("all")
    with Image.open(path) as image:
        assert image.size == (int(30 / 25.4 * 100),) * 2


def test_A1_resolve_alias():
    assert _resolve_alias("a", None, "new", "old") == "a"
    assert _resolve_alias(None, "b", "new", "old") == "b"
    assert _resolve_alias("a", "a", "new", "old") == "a"
    with pytest.raises(TypeError, match="aliases"):
        _resolve_alias("a", "b", "new", "old")
