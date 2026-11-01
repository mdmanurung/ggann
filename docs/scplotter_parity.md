# scplotter comparisons

This page maps scplotter's `CellDimPlot` and `FeatureStatPlot` interfaces to
ggann helpers and plotnine layers. The rendered examples use
`pbmc3k_processed`; scplotter's vignettes use a different Seurat dataset, so the
comparison concerns plot types and options rather than pixel equality.

[`examples/reproduce_scplotter.py`](https://github.com/mdmanurung/ggann/blob/main/examples/reproduce_scplotter.py)
generates the images.

```python
import scanpy as sc
import ggann as ag

adata = sc.datasets.pbmc3k_processed()
group = "louvain"
genes = ["CD3D", "NKG7", "CST3"]
```

## CellDimPlot

| scplotter | ggann |
|---|---|
| `CellDimPlot(group_by=, reduction="UMAP")` | `ag.plot_embedding(adata, "umap", color=group)` |
| `CellDimPlot(..., label=TRUE)` | `ag.plot_embedding(adata, "umap", color=group, label=True)` |
| `CellDimPlot(..., split_by="Phase")` | `ag.plot_embedding(adata, "umap", color=group, split_by="phase")` |
| `FeatureStatPlot(plot_type="dim", features=)` | `ag.plot_embedding(adata, "umap", color="CD3D")` |

| | |
|:---:|:---:|
| ![Basic embedding](images/scplotter/celldim_basic.png) | ![Embedding labels](images/scplotter/celldim_label.png) |
| ![Split embedding](images/scplotter/celldim_split.png) | ![Feature embedding](images/scplotter/celldim_feature.png) |

Density contours and binned counts use ordinary plotnine layers:

```python
from plotnine import aes, geom_bin2d, geom_density_2d, ggplot

density = ag.plot_embedding(adata, "umap", color=group) + geom_density_2d(
    color="black"
)

coords = ag.embedding_coords(adata, "umap")
x, y = coords.columns[:2]
hexbin = ggplot(coords, aes(x, y)) + geom_bin2d(bins=28)
```

| | |
|:---:|:---:|
| ![Density contours](images/scplotter/celldim_density.png) | ![Binned embedding](images/scplotter/celldim_hex.png) |

Features that require lineage fits, velocity estimates, neighbour graphs,
three-dimensional rendering, or embedded charts are outside ggann's data and
plotting model. `plot_embedding` has no `projection="3d"`, `edges`, `arrows` or
`add_outline`; see {doc}`limitations`.

## FeatureStatPlot

| scplotter option | ggann |
|---|---|
| Default violin | `ag.plot_violin(adata, genes, group)` |
| `add_point=TRUE` | `ag.plot_violin(..., add_points=True)` |
| `plot_type="box"` | `ag.plot_box(adata, genes, group)` |
| `plot_type="bar"` | `ag.plot_expression_bar(adata, genes, group)` |
| `stack=TRUE` | `ag.plot_stacked_violin(adata, genes, group)` |
| `comparisons=TRUE` | `ag.plot_violin(..., stats=True)` |
| `plot_type="heatmap"` | `ag.plot_matrixplot(adata, genes, group)` |
| `plot_type="dot"` | `ag.plot_dotplot(adata, genes, group)` |
| `plot_type="cor"` | `ag.plot_correlation(adata, group, genes=genes)` |
| `plot_type="ridge"` | `ag.plot_ridge(adata, genes, group)` |

| | |
|:---:|:---:|
| ![Violin](images/scplotter/feat_violin.png) | ![Violin with points](images/scplotter/feat_violin_points.png) |
| ![Box](images/scplotter/feat_box.png) | ![Expression bar](images/scplotter/feat_bar.png) |
| ![Stacked violin](images/scplotter/feat_stacked_violin.png) | ![Statistical annotation](images/scplotter/feat_stats.png) |
| ![Matrix plot](images/scplotter/feat_heatmap.png) | ![Dotplot](images/scplotter/feat_dot.png) |
| ![Correlation](images/scplotter/feat_cor.png) | ![Flipped boxplot](images/scplotter/feat_box_flip.png) |

The stacked violin fills each violin by group median expression, as Scanpy
does, and draws genes as rows with no y tick labels. `plot_violin`,
`plot_box` and `plot_stacked_violin` take `order=` to subset and order the
groups, and `log=True` to plot `log1p(expression)`. The dotplot and matrixplot
draw the first category at the top.
`stat_compare_means(comparisons=)` fails upstream, so do not use it for pairwise
brackets; see {doc}`limitations`.

Custom summaries and coordinate changes remain plotnine operations:

```python
import numpy as np
from plotnine import coord_flip, stat_summary

with_mean = ag.plot_violin(adata, ["CD3D"], group, add_box=False) + stat_summary(
    fun_y=np.mean,
    geom="point",
)
flipped = ag.plot_box(adata, ["CD3D"], group) + coord_flip()
```

Annotated grid heatmaps can use `plot_clustermap`; correlation pairs matrices
and lineage-aware feature plots are not implemented.
Use `swap_axes=True`, not `coord_flip()`, to swap the dotplot or matrixplot
axes, because a flip reverses the displayed category order. `coord_flip()` is
fine for the box plot above. The images on this page were regenerated with
`python examples/reproduce_scplotter.py`.
