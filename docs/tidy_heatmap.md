# tidyHeatmap-style heatmaps

`tidy_heatmap` is a Python port of the verbs of the R package
[tidyHeatmap](https://github.com/stemangiola/tidyHeatmap), drawn with
PyComplexHeatmap. Each verb returns a new immutable `TidyHeatmap`. Nothing is
drawn until `.draw()`, `.as_ComplexHeatmap()`, `.save()` or Jupyter display.
Install the backend with `pip install ggann[heatmap]`.

Like `plot_clustermap`, the result is a grid-backend object, not a plotnine
object, so `compose` does not accept it. Use `+` to place heatmaps side by side.

## Two inputs

A long table needs a `row`, a `column` and a `value` column. A
`DataFrameGroupBy` plays the role of R's `group_by()`. An `AnnData` needs
`genes`: rows are genes, columns are cells or, with `group_by`, group means.

Every other column of a long table is routed like tidyHeatmap does. A column
that is constant within each row entity annotates the rows. One that is constant
within each column entity annotates the columns. A column that is constant
within neither raises a `ValueError` that names it. For an `AnnData`, `var`
columns annotate genes and `obs` columns annotate cells (with `group_by`, only
`obs` columns that are constant within each group, plus `n_cells`).

## R to Python

| tidyHeatmap (R) | ggann (Python) | Support |
|---|---|---|
| `heatmap(.row, .column, .value, transform, scale, palette_value, palette_grouping)` | `tidy_heatmap(df, row, column, value, transform=, scale=, palette_value=, palette_grouping=)` | Full. `na=` is new: missing values raise unless you pass a policy. |
| `group_by(x) \|> heatmap(...)` | `tidy_heatmap(df.groupby("x"), ...)` | Full. One grouping per axis, as in R. |
| (AnnData input) | `tidy_heatmap(adata, genes=, group_by=, layer=, use_raw=)` | New in ggann. |
| `annotation_tile(.column, palette, size)` | `.annotation_tile(column, palette=, size=)` | Full. `size` is in millimetres. |
| `annotation_point` | `.annotation_point(column)` | Full. |
| `annotation_bar` | `.annotation_bar(column)` | Full. Bars start at 0. |
| `annotation_line` | `.annotation_line(column)` | Full. |
| `annotation_numeric(labels_format)` | `.annotation_numeric(column, labels_format=)` | Partial. PyComplexHeatmap has no `anno_numeric`; this is a bar annotation with the value printed in each cell. |
| `annotation_group(..., group_label_fontsize, show_group_name, group_strip_height)` | `.annotation_group("x", ...)` | Partial. `show_group_name=False` blanks the strip title instead of removing it. |
| `layer_point(..., .size)` | `.layer_point("padj < 0.05", size=)` | Full. The filter is a pandas `query` string or a callable returning a mask. |
| `layer_square`, `layer_diamond`, `layer_arrow_up`, `layer_arrow_down`, `layer_star`, `layer_asterisk` | same names | Full. Star is R symbol 8 (8 spokes), asterisk is symbol 42 (6 spokes). |
| `layer_text(..., .value, .size, .color)` | `.layer_text("padj < 0.05", label="stars")` | Full. |
| `split_rows(k)`, `split_columns(k)` | same names | Full. Cuts the dendrogram with the heatmap's own method and metric. |
| `get_heatmap_data()` | `.get_heatmap_data()` | Full for the matrix and orders. R's dendrogram objects are not returned. |
| `as_ComplexHeatmap()` | `.as_ComplexHeatmap()` | Returns a `PyComplexHeatmap.ClusterMapPlotter`. |
| `a + b` | `a + b` | Partial. Uses `PyComplexHeatmap.composite`; see below. |
| `save_pdf()` | `.save(path, width=, height=)` | Delegates to `save_publication`. |

## Differences to know

- Clustering defaults match ComplexHeatmap (`complete`, Euclidean), not
  PyComplexHeatmap's `ward`. Override with `row_cluster_method=`.
- Value colours interpolate `palette_value` in RGB; R interpolates in LAB.
- Scaling is `scale_robust`: a z-score with the sample standard deviation, where
  a standard deviation of 0 gives zeros instead of NaN. `transform` runs first.
- Duplicated (row, column) pairs raise even if their values agree.
- With `+`, all heatmaps need the same rows, and the first one fixes their order.
  PyComplexHeatmap does not align the heatmap bodies when annotation heights
  differ, so give each heatmap the same number of column annotations.
- Extra keyword arguments go to `ClusterMapPlotter`, for example
  `row_cluster=False` or `show_rownames=False`.
- Size the figure at export: `.save(...)` draws at the export size. Long labels
  and legends can still overflow a small canvas; `save_publication` warns.

## Example: a long differential-expression table

```python
import numpy as np
import pandas as pd
import ggann as ag

rng = np.random.default_rng(0)
lineage = {"T": "lymphoid", "B": "lymphoid", "NK": "lymphoid", "Mono": "myeloid", "DC": "myeloid"}
df = pd.DataFrame(
    [
        {
            "feature": f"G{i:02d}",
            "cell_type": cell_type,
            "mean": rng.gamma(2.0, 1.0),
            "padj": rng.random() ** 2,
            "lineage": lineage[cell_type],
            "n_cells": 20 * (k + 1),
        }
        for i in range(14)
        for k, cell_type in enumerate(lineage)
    ]
)

heatmap = (
    ag.tidy_heatmap(df, row="feature", column="cell_type", value="mean", scale="row")
    .annotation_tile("lineage")
    .annotation_bar("n_cells")
    .split_rows(3)
    .layer_point("padj < 0.05")
)
plotter = heatmap.draw(figsize=(7, 6))
heatmap.get_heatmap_data().matrix.head()
```

## Example: AnnData with `group_by`

```python
import scanpy as sc
import ggann as ag

adata = sc.datasets.pbmc68k_reduced()
heatmap = (
    ag.tidy_heatmap(
        adata,
        genes=["CD3D", "NKG7", "CST3", "MS4A1", "GNLY", "LYZ"],
        group_by="bulk_labels",
        scale="row",
        show_colnames=False,  # the tile strip names the groups
        plot_legend=False,  # legends overflow a 183 mm canvas; draw them separately
    )
    .annotation_tile("bulk_labels")
    .annotation_bar("n_cells")
    .layer_point("mean_expression > 1", size=2)
)
heatmap.save("pbmc_heatmap", width="double-column", height=90, formats=("png",))
```

Grouping with `df.groupby("lineage")` splits the columns and adds a group strip:

```python
grouped = ag.tidy_heatmap(
    df.groupby("lineage"), row="feature", column="cell_type", value="mean", scale="row"
)
```
