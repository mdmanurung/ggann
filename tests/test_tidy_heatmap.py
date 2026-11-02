"""tidyHeatmap-style interface: verbs, routing, scaling, validation and export."""

from __future__ import annotations

import subprocess
import sys
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

import ggann as ag
from ggann.tidy_heatmap import TidyHeatmapSet

pytest.importorskip("PyComplexHeatmap")

LINEAGE = {"T": "lymphoid", "B": "lymphoid", "NK": "lymphoid", "Mono": "myeloid", "DC": "myeloid"}


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    plt.close("all")


@pytest.fixture
def long_df():
    """Six features x five cell types; padj/lineage/n_cells/gene_len route to an axis."""
    rng = np.random.default_rng(0)
    rows = [
        {
            "feature": f"G{i}",
            "cell_type": ct,
            "mean": float(rng.gamma(2.0, 1.0)),
            "padj": float(rng.random()),
            "lineage": LINEAGE[ct],
            "n_cells": 10.0 * (k + 1),
            "module": f"m{i % 2}",
            "gene_len": float(1 + i),
        }
        for i in range(6)
        for k, ct in enumerate(LINEAGE)
    ]
    return pd.DataFrame(rows)


def _hm(df, **kwargs):
    return ag.tidy_heatmap(df, row="feature", column="cell_type", value="mean", **kwargs)


def _matrix(hm):
    return hm._spec.matrix


# ------------------------------------------------------------------ package contract


def test_import_does_not_need_pycomplexheatmap():
    code = "import sys; import ggann; assert 'PyComplexHeatmap' not in sys.modules"
    subprocess.run([sys.executable, "-c", code], check=True)
    assert ag.tidy_heatmap is not None and ag.TidyHeatmap is not None
    assert {"tidy_heatmap", "TidyHeatmap"} <= set(ag.__all__)


# ----------------------------------------------------------------- scaling numerics


def _zscore(frame, axis):
    mean, sd = frame.mean(axis=axis), frame.std(axis=axis, ddof=1)
    sd = sd.where(sd != 0, 1.0)
    return frame.sub(mean, axis=1 - axis).div(sd, axis=1 - axis)


def test_scale_row_matches_manual_zscore_with_a_constant_row(long_df):
    long_df.loc[long_df.feature == "G0", "mean"] = 3.0  # sd == 0
    wide = long_df.pivot(index="feature", columns="cell_type", values="mean")
    got = _matrix(_hm(long_df, scale="row"))
    assert not got.isna().any().any()
    assert (got.loc["G0"] == 0).all()
    pd.testing.assert_frame_equal(got, _zscore(wide, 1), check_names=False)


def test_scale_column_and_both_match_manual_zscore(long_df):
    wide = long_df.pivot(index="feature", columns="cell_type", values="mean")
    pd.testing.assert_frame_equal(
        _matrix(_hm(long_df, scale="column")), _zscore(wide, 0), check_names=False
    )
    # both = rows first, then columns, as tidyHeatmap
    pd.testing.assert_frame_equal(
        _matrix(_hm(long_df, scale="both")), _zscore(_zscore(wide, 1), 0), check_names=False
    )


def test_scale_none_keeps_values(long_df):
    wide = long_df.pivot(index="feature", columns="cell_type", values="mean")
    pd.testing.assert_frame_equal(_matrix(_hm(long_df)), wide, check_names=False)


def test_scale_of_a_single_value_is_zero_not_nan():
    df = pd.DataFrame({"r": ["a", "b"], "c": ["x", "x"], "v": [1.0, 5.0]})
    got = _matrix(ag.tidy_heatmap(df, row="r", column="c", value="v", scale="row"))
    assert (got.to_numpy() == 0).all()


def test_transform_runs_before_scaling(long_df):
    wide = long_df.pivot(index="feature", columns="cell_type", values="mean")
    got = _matrix(_hm(long_df, transform=np.log1p, scale="row"))
    pd.testing.assert_frame_equal(got, _zscore(np.log1p(wide), 1), check_names=False)
    # scaling first and transforming later would give a different matrix
    assert not np.allclose(got, np.log1p(_zscore(wide, 1).clip(lower=0)))


def test_invalid_scale_and_transform(long_df):
    with pytest.raises(ValueError, match="scale must be one of"):
        _hm(long_df, scale="rows")
    with pytest.raises(TypeError, match="callable"):
        _hm(long_df, transform="log1p")
    with pytest.raises(ValueError, match="negative infinite"):
        _hm(long_df.assign(mean=0.0), transform=np.log)
    with pytest.raises(ValueError, match="introduced NaN"):
        _hm(long_df.assign(mean=-1.0), transform=np.sqrt)


# ----------------------------------------------------------------------- validation


def test_duplicated_pairs_raise(long_df):
    with pytest.raises(ValueError, match="duplicated feature/cell_type pairs"):
        _hm(pd.concat([long_df, long_df.head(1)]))


def test_missing_values_raise_unless_a_policy_is_given(long_df):
    holey = long_df.copy()
    holey.loc[0, "mean"] = np.nan
    with pytest.raises(ValueError, match="1 missing cells"):
        _hm(holey)
    with pytest.warns(UserWarning, match="dropped 1 of 6 feature"):
        dropped = _hm(holey, na="drop")
    assert dropped._spec.matrix.shape == (5, 5)
    assert _matrix(_hm(holey, na=0.0)).loc["G0", "T"] == 0.0
    absent_pair = long_df.iloc[1:]
    with pytest.raises(ValueError, match="missing cells"):
        _hm(absent_pair)


def test_bad_columns_and_values(long_df):
    with pytest.raises(KeyError, match="nope"):
        ag.tidy_heatmap(long_df, row="feature", column="cell_type", value="nope")
    with pytest.raises(ValueError, match="must be numeric"):
        ag.tidy_heatmap(long_df, row="feature", column="cell_type", value="lineage")
    with pytest.raises(ValueError, match="palette_value"):
        _hm(long_df, palette_value=["#000000", "#ffffff"])
    with pytest.raises(ValueError, match="use palette_value"):
        _hm(long_df, cmap="viridis")
    with pytest.raises(ValueError, match="only for AnnData"):
        _hm(long_df, genes=["G0"])
    with pytest.raises(ValueError, match="genes"):
        ag.tidy_heatmap(pd.DataFrame(), row="a", column="b", value="c", genes=["x"])


# -------------------------------------------------------------------------- routing


def test_annotations_route_to_rows_or_columns_by_constancy(long_df):
    hm = (
        _hm(long_df).annotation_tile("lineage").annotation_tile("module").annotation_bar("gene_len")
    )
    axes = {a.column: a.axis for a in hm._spec.annotations}
    assert axes == {"lineage": "column", "module": "row", "gene_len": "row"}
    plotter = hm.draw()
    assert [a.df.columns[0] for a in plotter.top_annotation.annotations] == ["lineage"]
    assert [a.df.columns[0] for a in plotter.left_annotation.annotations] == ["module", "gene_len"]


def test_a_column_constant_within_neither_axis_raises_and_names_it(long_df):
    with pytest.raises(ValueError, match="'padj'"):
        _hm(long_df).annotation_tile("padj")
    with pytest.raises(KeyError, match="'absent'"):
        _hm(long_df).annotation_tile("absent")
    with pytest.raises(ValueError, match="already an annotation"):
        _hm(long_df).annotation_tile("lineage").annotation_tile("lineage")
    with pytest.raises(ValueError, match="needs a numeric column"):
        _hm(long_df).annotation_bar("lineage")


def test_the_group_by_column_itself_is_a_valid_annotation(long_df):
    hm = _hm(long_df).annotation_tile("cell_type")
    assert hm._spec.annotations[0].axis == "column"
    hm.draw()


# --------------------------------------------------------------------- every verb


@pytest.mark.parametrize(
    "verb,column",
    [
        ("annotation_tile", "lineage"),
        ("annotation_tile", "n_cells"),
        ("annotation_point", "n_cells"),
        ("annotation_bar", "n_cells"),
        ("annotation_line", "n_cells"),
        ("annotation_numeric", "n_cells"),
        ("annotation_bar", "gene_len"),
    ],
)
def test_annotation_verbs_draw_on_a_long_table(long_df, verb, column):
    hm = getattr(_hm(long_df), verb)(column)
    plotter = hm.draw()
    holder = plotter.left_annotation if column == "gene_len" else plotter.top_annotation
    assert holder is not None and holder.annotations[0].df.columns[0] == column


def test_annotation_numeric_prints_the_values(long_df):
    plotter = _hm(long_df).annotation_numeric("n_cells", labels_format="{:.0f}").draw()
    ax = plotter.top_annotation.annotations[0].ax
    assert sorted(t.get_text() for t in ax.texts) == ["10", "20", "30", "40", "50"]


def test_annotation_palette_and_size(long_df):
    palette = {"lymphoid": "#ff0000", "myeloid": "#0000ff"}
    plotter = _hm(long_df).annotation_tile("lineage", palette=palette, size=6).draw()
    anno = plotter.top_annotation.annotations[0]
    assert {k: matplotlib_hex(v) for k, v in anno.color_dict.items()} == palette
    assert anno.height == 6


def matplotlib_hex(color):
    from matplotlib.colors import to_hex

    return to_hex(color)


def test_default_discrete_palette_follows_tidyheatmap_set1(long_df):
    from matplotlib import colormaps
    from matplotlib.colors import to_hex

    plotter = _hm(long_df).annotation_tile("lineage").draw()
    colours = list(plotter.top_annotation.annotations[0].color_dict.values())
    set1 = [to_hex(c) for c in colormaps["Set1"].colors]
    # character columns interpolate the whole palette: first and last of Set1 for two levels
    assert [to_hex(c) for c in colours] == [set1[0], set1[-1]]


def _scatter_count(plotter):
    return sum(len(ax.collections) for ax in plotter.heatmap_axes.ravel()) - len(
        plotter.heatmap_axes.ravel()
    )


@pytest.mark.parametrize(
    "verb",
    [
        "layer_point",
        "layer_square",
        "layer_diamond",
        "layer_star",
        "layer_asterisk",
        "layer_arrow_up",
        "layer_arrow_down",
    ],
)
def test_layer_symbol_verbs_mark_the_selected_cells(long_df, verb):
    n = int((long_df.padj < 0.5).sum())
    for where in ("padj < 0.5", lambda d: d["padj"] < 0.5):
        plotter = getattr(_hm(long_df), verb)(where, size=2, color="red").draw()
        assert _scatter_count(plotter) == n


def test_layer_marks_land_on_the_cell_they_describe(long_df):
    df = long_df.assign(hit=(long_df.feature == "G2") & (long_df.cell_type == "NK"))
    plotter = _hm(df).layer_point("hit").draw()
    (ax,) = [a for a in plotter.heatmap_axes.ravel() if len(a.collections) > 1]
    x, y = ax.collections[-1].get_offsets()[0]
    rows = [r for part in plotter.row_order for r in part]
    cols = [c for part in plotter.col_order for c in part]
    assert (rows[int(y)], cols[int(x)]) == ("G2", "NK")


def test_layer_text_prints_the_label_column(long_df):
    df = long_df.assign(stars=np.where(long_df.padj < 0.3, "*", ""))
    n = int((long_df.padj < 0.3).sum())
    plotter = _hm(df).layer_text("padj < 0.3", label="stars", size=7, color="white").draw()
    texts = [t for ax in plotter.heatmap_axes.ravel() for t in ax.texts]
    assert len(texts) == n and {t.get_text() for t in texts} == {"*"}
    with pytest.raises(KeyError, match="'nope'"):
        _hm(df).layer_text(label="nope")
    with pytest.raises(TypeError):
        _hm(df).layer_point(5)  # type: ignore[arg-type]


def test_layer_size_and_color_can_come_from_columns(long_df):
    df = long_df.assign(sz=2.0, col="#ff0000")
    hm = _hm(df).layer_point("padj < 0.5", size="sz", color="col")
    assert set(hm._spec.layers[0].cells["color"]) == {"#ff0000"}
    hm.draw()


@pytest.mark.parametrize("axis", ["rows", "columns"])
def test_split_verbs_cut_the_dendrogram(long_df, axis):
    hm = getattr(_hm(long_df), f"split_{axis}")(2)
    plotter = hm.draw()
    order = plotter.row_order if axis == "rows" else plotter.col_order
    assert len(order) == 2
    with pytest.raises(ValueError, match="between 2 and"):
        getattr(_hm(long_df), f"split_{axis}")(99)


def test_split_rows_matches_a_scipy_cut_of_the_same_tree(long_df):
    from scipy.cluster.hierarchy import fcluster, linkage

    hm = _hm(long_df, scale="row").split_rows(3)
    plotter = hm.draw()
    m = _matrix(hm)
    labels = fcluster(linkage(m.to_numpy(), "complete", "euclidean"), 3, "maxclust")
    expected = {frozenset(m.index[labels == k]) for k in set(labels)}
    assert {frozenset(part) for part in plotter.row_order} == expected


def test_split_conflicts_are_reported(long_df):
    with pytest.raises(ValueError, match="clustering is switched off"):
        _hm(long_df, row_cluster=False).split_rows(2)
    grouped = ag.tidy_heatmap(
        long_df.groupby("module"), row="feature", column="cell_type", value="mean"
    )
    with pytest.raises(ValueError, match="already split"):
        grouped.split_rows(2)


# ---------------------------------------------------------------------- group_by()


def test_groupby_routes_row_and_column_groups(long_df):
    hm = ag.tidy_heatmap(
        long_df.groupby(["lineage", "module"]), row="feature", column="cell_type", value="mean"
    )
    assert (hm._spec.row_group, hm._spec.col_group) == ("module", "lineage")
    plotter = hm.draw()
    assert len(plotter.row_order) == 2 and len(plotter.col_order) == 2
    # slices follow sorted group labels (cluster_row_slices = FALSE)
    assert set(plotter.row_order[0]) == {"G0", "G2", "G4"}
    assert set(plotter.col_order[0]) == {"B", "NK", "T"}
    assert [a.df.columns[0] for a in plotter.top_annotation.annotations] == ["lineage"]


def test_groupby_strips_use_palette_grouping(long_df):
    hm = ag.tidy_heatmap(
        long_df.groupby("lineage"),
        row="feature",
        column="cell_type",
        value="mean",
        palette_grouping=[["#ff0000", "#00ff00"]],
    )
    anno = hm.draw().top_annotation.annotations[0]
    assert {matplotlib_hex(c) for c in anno.color_dict.values()} == {"#ff0000", "#00ff00"}


def test_annotation_group_matches_groupby(long_df):
    by_verb = _hm(long_df).annotation_group("lineage", group_label_fontsize=6)
    by_frame = ag.tidy_heatmap(
        long_df.groupby("lineage"), row="feature", column="cell_type", value="mean"
    )
    assert by_verb._spec.col_group == by_frame._spec.col_group == "lineage"
    assert by_verb._spec.group_label_fontsize == 6
    by_verb.draw()


def test_groupby_errors(long_df):
    with pytest.raises(ValueError, match="'padj'"):
        ag.tidy_heatmap(long_df.groupby("padj"), row="feature", column="cell_type", value="mean")
    with pytest.raises(ValueError, match="one grouping per dimension"):
        ag.tidy_heatmap(
            long_df.assign(m2="z").groupby(["module", "gene_len"]),
            row="feature",
            column="cell_type",
            value="mean",
        )


# --------------------------------------------------------- heatmap data and drawing


def test_get_heatmap_data_matches_the_drawn_plotter(long_df):
    hm = _hm(long_df, scale="row").split_rows(2).split_columns(2)
    data = hm.get_heatmap_data()
    plotter = hm.as_ComplexHeatmap()
    assert data.row_order == [r for part in plotter.row_order for r in part]
    assert data.column_order == [c for part in plotter.col_order for c in part]
    pd.testing.assert_frame_equal(
        data.matrix, plotter.data2d.loc[data.row_order, data.column_order]
    )
    assert list(data.matrix.index) == data.row_order
    assert plt.get_fignums() == [plotter.ax.figure.number]  # get_heatmap_data left no figure


def test_unclustered_heatmap_keeps_sorted_order(long_df):
    data = _hm(long_df, row_cluster=False, col_cluster=False).get_heatmap_data()
    assert data.row_order == sorted(long_df.feature.unique())
    assert data.column_order == sorted(long_df.cell_type.unique())


def test_categorical_keys_keep_their_level_order(long_df):
    order = ["NK", "T", "DC", "B", "Mono"]
    df = long_df.assign(cell_type=pd.Categorical(long_df.cell_type, categories=order))
    data = _hm(df, col_cluster=False).get_heatmap_data()
    assert data.column_order == order


def test_value_colours_follow_palette_value_and_constant_matrix():
    df = pd.DataFrame({"r": list("abab"), "c": list("xxyy"), "v": [2.0] * 4})
    plotter = ag.tidy_heatmap(df, row="r", column="c", value="v").draw()
    assert (plotter.kwargs["vmin"], plotter.kwargs["vmax"]) == (1.0, 3.0)  # tidyHeatmap issue 40
    assert plotter.cmap.N == 256 and matplotlib_hex(plotter.cmap(0.0)) == "#440154"


def test_clustering_needs_two_entries():
    df = pd.DataFrame({"r": ["a"], "c": ["x"], "v": [1.0]})
    with pytest.raises(ValueError, match="row_cluster=False"):
        ag.tidy_heatmap(df, row="r", column="c", value="v").draw()
    ag.tidy_heatmap(df, row="r", column="c", value="v", row_cluster=False, col_cluster=False).draw()


# ----------------------------------------------------------- immutability and input


def test_verbs_return_new_builders_and_the_input_is_untouched(long_df):
    before = long_df.copy()
    base = _hm(long_df)
    chained = base.annotation_tile("lineage").layer_point("padj < 0.1").split_rows(2)
    assert chained is not base
    assert base._spec.annotations == () and base._spec.layers == ()
    assert base._spec.row_split is None and chained._spec.row_split == 2
    with pytest.raises(AttributeError, match="immutable"):
        base.anything = 1
    chained.draw()
    pd.testing.assert_frame_equal(long_df, before)
    long_df.loc[0, "mean"] = 999.0  # later edits to the input do not leak into the builder
    assert _matrix(base).max().max() < 999.0


def test_repr_png_returns_a_png_and_leaks_no_figure(long_df):
    png = _hm(long_df)._repr_png_()
    assert png[:8] == b"\x89PNG\r\n\x1a\n" and plt.get_fignums() == []


# ------------------------------------------------------------------------ AnnData


GENES = ["CD3D", "NKG7", "CST3", "MS4A1"]


@pytest.fixture
def genes(adata):
    return [g for g in GENES if g in adata.raw.var_names]


def test_anndata_group_means_match_group_means(adata, genes, group_key):
    from ggann._aggregate import group_means

    hm = ag.tidy_heatmap(adata, genes=genes, group_by=group_key)
    expected = group_means(adata, genes, group_key).T.loc[genes]
    got = _matrix(hm)
    assert list(got.index) == genes
    np.testing.assert_allclose(got.to_numpy(), expected.to_numpy())
    assert hm._spec.value == "mean_expression"


def test_anndata_scale_and_groupby_alias(adata, genes, group_key):
    a = ag.tidy_heatmap(adata, genes=genes, group_by=group_key, scale="row")
    b = ag.tidy_heatmap(adata, genes=genes, groupby=group_key, scale="row")
    pd.testing.assert_frame_equal(_matrix(a), _matrix(b))
    raw = ag.tidy_heatmap(adata, genes=genes, group_by=group_key)
    pd.testing.assert_frame_equal(_matrix(a), _zscore(_matrix(raw), 1), check_names=False)
    with pytest.raises(TypeError, match="aliases"):
        ag.tidy_heatmap(adata, genes=genes, group_by=group_key, groupby="louvain")


def test_anndata_cells_as_columns_with_obs_and_var_annotations(adata, genes):
    hm = (
        ag.tidy_heatmap(adata, genes=genes, use_raw=False, scale="row")
        .annotation_tile("bulk_labels")
        .annotation_bar("n_counts")
        .annotation_tile("highly_variable")
        .annotation_bar("means")
    )
    assert {a.column: a.axis for a in hm._spec.annotations} == {
        "bulk_labels": "column",
        "n_counts": "column",
        "highly_variable": "row",
        "means": "row",
    }
    plotter = hm.draw()
    assert plotter.data2d.shape == (len(genes), adata.n_obs)
    assert plotter.left_annotation is not None
    # stored obs colours are reused
    stored = dict(zip(adata.obs.bulk_labels.cat.categories, adata.uns["bulk_labels_colors"]))
    tile = plotter.top_annotation.annotations[0].color_dict
    assert {k: matplotlib_hex(v) for k, v in tile.items()} == stored


def test_anndata_group_by_annotations_and_n_cells(adata, genes, group_key):
    adata = adata.copy()
    adata.obs["lineage"] = adata.obs[group_key].astype(str).str[:1]  # constant within a group
    hm = (
        ag.tidy_heatmap(adata, genes=genes, group_by=group_key)
        .annotation_tile("lineage")
        .annotation_bar("n_cells")
        .annotation_tile(group_key)
    )
    assert all(a.axis == "column" for a in hm._spec.annotations)
    counts = adata.obs[group_key].value_counts()
    attrs = hm._spec.col_attrs
    assert (attrs["n_cells"] == counts.reindex(attrs.index)).all()
    hm.draw()
    # an obs column that varies within a group cannot be routed
    with pytest.raises(ValueError, match="'n_counts'"):
        ag.tidy_heatmap(adata, genes=genes, group_by=group_key).annotation_tile("n_counts")


def test_anndata_layers_splits_and_group_strips(adata, genes, group_key):
    hm = (
        ag.tidy_heatmap(adata, genes=genes, group_by=group_key, scale="row")
        .layer_diamond("mean_expression > 1", color="white")
        .layer_text("mean_expression > 1", label="gene", size=5)
        .split_rows(2)
        .split_columns(3)
    )
    plotter = hm.draw()
    assert len(plotter.row_order) == 2 and len(plotter.col_order) == 3
    n = int((_matrix(hm) > _matrix(hm).min().min()).sum().sum())  # layers drew something
    assert n > 0 and _scatter_count(plotter) > 0
    stripped = ag.tidy_heatmap(adata, genes=genes, group_by=group_key).annotation_group(group_key)
    assert stripped._spec.col_group == group_key
    stripped.draw()


def test_anndata_input_is_not_modified_and_validates(adata, genes, group_key):
    obs_before, x_before = adata.obs.copy(), adata.X.copy()
    (
        ag.tidy_heatmap(adata, genes=genes, group_by=group_key, scale="both")
        .annotation_tile(group_key)
        .layer_point("mean_expression > 0")
        .draw()
    )
    pd.testing.assert_frame_equal(adata.obs, obs_before)
    assert (adata.X == x_before).all()
    with pytest.raises(ValueError, match="genes"):
        ag.tidy_heatmap(adata)
    with pytest.raises(ValueError, match="long tables"):
        ag.tidy_heatmap(adata, "gene", genes=genes)
    with pytest.raises(KeyError):
        ag.tidy_heatmap(adata, genes=["NOT_A_GENE"])
    with pytest.raises(KeyError):
        ag.tidy_heatmap(adata, genes=genes, group_by="nope")


def test_anndata_drops_cells_with_a_missing_group(adata, genes, group_key):
    adata = adata.copy()
    adata.obs[group_key] = adata.obs[group_key].astype(object)
    adata.obs.iloc[:3, adata.obs.columns.get_loc(group_key)] = np.nan
    with pytest.warns(UserWarning, match="dropped 3 of"):
        ag.tidy_heatmap(adata, genes=genes, group_by=group_key)


# ------------------------------------------------------------------ + and saving


def test_plus_places_heatmaps_side_by_side(long_df):
    a = _hm(long_df, scale="row").annotation_tile("lineage")
    b = (
        _hm(long_df.assign(mean=long_df["mean"] * 2))
        .annotation_bar("n_cells")
        .layer_point("padj < 0.2")
    )
    both = a + b + a
    assert isinstance(both, TidyHeatmapSet)
    assert len(both.heatmaps) == 3
    fig = both.draw()
    assert sum(len(ax.collections) > 0 for ax in fig.axes) >= 3
    assert both._repr_png_()[:4] == b"\x89PNG"
    other_rows = long_df[long_df.feature != "G0"]
    with pytest.raises(ValueError, match="same rows"):
        _ = a + _hm(other_rows)


def test_save_goes_through_save_publication(long_df, tmp_path):
    hm = _hm(long_df, scale="row").annotation_tile("lineage").layer_point("padj < 0.2")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)  # PyComplexHeatmap overflows small canvases
        paths = hm.save(tmp_path / "fig", width=89, height=70, formats=("png", "svg"), dpi=150)
    assert [p.suffix for p in paths] == [".svg", ".png"]
    image = plt.imread(tmp_path / "fig.png")
    assert image.shape[:2] == (int(70 / 25.4 * 150), int(89 / 25.4 * 150))
    with pytest.raises(ValueError, match="height must be a positive number"):
        hm.save(tmp_path / "bad", height=-1)


def test_save_of_a_composition(long_df, tmp_path):
    both = _hm(long_df) + _hm(long_df)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        (path,) = both.save(tmp_path / "two.png", width=183, height=80, dpi=100)
    assert path.exists()
