# ggann

`ggann` turns `AnnData` into composable plotnine graphics and concise
single-cell plotting helpers.

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

Run five offline workflows from migration through final figure assembly.
:::

:::{grid-item-card} Performance
:link: performance
:link-type: doc

Reproduce matched Scanpy timings, memory measurements, and release gates.
:::

::::

```{toctree}
:caption: Get started
:maxdepth: 1
:hidden:

installation
quickstart
concepts
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
:caption: Release evidence
:maxdepth: 1
:hidden:

performance
migration
```

```{toctree}
:caption: Reference
:maxdepth: 1
:hidden:

api
```
