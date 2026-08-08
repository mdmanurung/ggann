# ggann

`ggann` prepares `AnnData` fields for plotnine and provides plotting helpers for
common single-cell figures. `gganndata(...)` and the plotnine-native helpers
return regular [`plotnine.ggplot`](https://plotnine.org/) objects.

```python
import scanpy as sc
import ggann as ag
from ggann import aes, gganndata
from plotnine import geom_point

adata = sc.datasets.pbmc68k_reduced()
markers = ["CD3D", "NKG7", "CST3"]

gganndata(
    adata,
    aes("UMAP_1", "UMAP_2", color="bulk_labels"),
) + geom_point(size=1.5)

ag.plot_dotplot(adata, markers, group_by="bulk_labels")
```

`plot_clustermap` and `plot_upset` use grid-based backends and return their
backend objects rather than `ggplot` objects.

## Start here

::::{grid} 1 1 2 2
:gutter: 2

:::{grid-item-card} Installation
:link: installation
:link-type: doc

Install the package and optional plotting backends.
:::

:::{grid-item-card} Quickstart
:link: quickstart
:link-type: doc

Build plots with the grammar and helper interfaces.
:::

:::{grid-item-card} Gallery
:link: gallery
:link-type: doc

View figures generated from `pbmc68k_reduced`.
:::

:::{grid-item-card} API reference
:link: api
:link-type: doc

Check signatures, parameters, and return types.
:::

::::

```{toctree}
:caption: Get Started
:maxdepth: 1
:hidden:

installation
quickstart
gallery
comparisons
scplotter_parity
scanpy_parity
stats_pseudobulk
```

```{toctree}
:caption: API Reference
:maxdepth: 1
:hidden:

api
```
