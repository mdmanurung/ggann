# ggann

`ggann` turns `AnnData` into composable plotnine graphics and concise
single-cell plotting helpers, then carries the same scientific content into an
exact-size editable publication figure.

```python
import ggann as ag

plot = ag.plot_embedding(adata, "umap", color="cell_type")
```

For custom figures, `gganndata` returns an ordinary `plotnine.ggplot`:

```python
from ggann import aes, gganndata
from plotnine import geom_point

plot = gganndata(
    adata,
    aes("UMAP_1", "UMAP_2", color="cell_type"),
) + geom_point()
```

All plotnine-native helpers remain composable. `plot_clustermap` and
`plot_upset` are the two documented grid-backend exceptions.

## Start here

::::{grid} 1 1 2 2
:gutter: 2

:::{grid-item-card} Quickstart
:link: quickstart
:link-type: doc

Create a helper plot, then build the same data flow with the grammar.
:::

:::{grid-item-card} API concepts
:link: concepts
:link-type: doc

Understand source resolution, grouping, ownership, and return types.
:::

:::{grid-item-card} Executable vignettes
:link: vignettes/index
:link-type: doc

Run six offline workflows, including two real PBMC analyses.
:::

:::{grid-item-card} Publication figures
:link: publication
:link-type: doc

Coordinate final-size styles, palettes, layouts, rasterization, and export.
:::

:::{grid-item-card} Performance
:link: performance
:link-type: doc

Reproduce matched Scanpy timings, memory measurements, and acceptance criteria.
:::

::::

```{toctree}
:caption: Get started
:maxdepth: 1
:hidden:

installation
quickstart
concepts
publication
gallery
comparisons
```

```{toctree}
:caption: Workflows
:maxdepth: 2
:hidden:

vignettes/index
scanpy_parity
scplotter_parity
stats_pseudobulk
```

```{toctree}
:caption: Adoption and validation
:maxdepth: 1
:hidden:

migration
performance
```

```{toctree}
:caption: Reference
:maxdepth: 1
:hidden:

api
```
