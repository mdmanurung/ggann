# Gallery

[`examples/gallery.py`](https://github.com/mdmanurung/ggann/blob/main/examples/gallery.py)
generates these figures from `pbmc68k_reduced`.

## Embeddings

::::{grid} 1 1 2 2

:::{grid-item}
![UMAP by cluster](images/umap_clusters.png)
:::

:::{grid-item}
![Labelled UMAP](images/umap_labelled.png)
:::

::::

`plot_embedding` colours cells by a categorical or numeric field. Set
`label=True` to add repelled category labels.

![Multi-gene grid](images/features_grid.png)

`plot_features` uses one panel per feature and a shared colour scale.

## Gene-weighted density

![Density](images/density.png)

`plot_density` uses pyNebulosa for weighted kernel density. Set `joint=True` to
add a joint panel.

## Markers and expression

::::{grid} 1 1 2 2

:::{grid-item}
![Dotplot](images/dotplot.png)
:::

:::{grid-item}
![Stacked violin](images/stacked_violin.png)
:::

:::{grid-item}
![Box](images/box.png)
:::

:::{grid-item}
![Expression bar](images/expression_bar.png)
:::

::::

The panels show `plot_dotplot`, `plot_stacked_violin`, `plot_box`, and
`plot_expression_bar`.

## Differential expression

::::{grid} 1 1 2 2

:::{grid-item}
![DE dotplot](images/de_dotplot.png)
:::

:::{grid-item}
![Volcano](images/volcano.png)
:::

::::

The panels show `plot_rank_genes_dotplot` and `plot_volcano`.

## Composition, correlation, and sets

::::{grid} 1 1 2 2

:::{grid-item}
![Composition](images/proportions.png)
:::

:::{grid-item}
![Correlation](images/correlation.png)
:::

:::{grid-item}
![UpSet](images/upset.png)
:::

:::{grid-item}
![Expression line](images/expression_line.png)
:::

::::

The panels show `plot_proportions`, `plot_correlation`, `plot_upset`, and
`plot_expression_line`.

## Quality control

![QC violins](images/qc_violin.png)

`plot_qc_violin` displays one facet per QC metric.
