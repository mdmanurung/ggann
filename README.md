# ggann

`ggann` provides a plotnine interface for `AnnData`. Its grammar entry point,
`gganndata`, returns a regular [`plotnine.ggplot`](https://plotnine.org/) object,
so plotnine layers, scales, facets, and themes compose with it.

[Documentation](https://mdmanurung.github.io/ggann/) ·
[Quickstart](https://mdmanurung.github.io/ggann/quickstart.html) ·
[API reference](https://mdmanurung.github.io/ggann/api.html) ·
[Gallery](https://mdmanurung.github.io/ggann/gallery.html)

## Quickstart

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

ag.plot_embedding(adata, basis="umap", color="CD3D")
ag.plot_dotplot(adata, markers, group_by="bulk_labels")
```

Plotnine-native helpers return `ggplot` objects. `plot_clustermap` and
`plot_upset` use grid-based backends and return their backend objects instead.

## Data resolution

`gganndata` resolves bare aesthetic names in this order:

| Order | Source | Example |
|---|---|---|
| 1 | `adata.obs` | `"bulk_labels"` |
| 2 | expression matrix | `"CD3D"` |
| 3 | `adata.obsm` coordinate | `"UMAP_1"` |

An observation column takes precedence when it has the same name as a gene.
Use a typed accessor or source prefix when the source should be explicit:

```python
from ggann import gene, obs, obsm

gganndata(
    adata,
    aes(
        x=obsm("umap", 0),
        y=obsm("umap", 1),
        color=gene("CD3D", use_raw=True),
        shape=obs("bulk_labels"),
    ),
) + geom_point()
```

Equivalent string prefixes are available:

| Prefix | Source | Example |
|---|---|---|
| `obs:` | `adata.obs` | `"obs:phase"` |
| `gene:` | plot-wide expression source | `"gene:CD3D"` |
| `gene:...@<layer>` | named layer | `"gene:CD3D@counts"` |
| `gene:...@raw` or `@X` | `adata.raw` or `adata.X` | `"gene:CD3D@raw"` |
| `obsm:<basis>[i]` | embedding coordinate, zero-based | `"obsm:umap[0]"` |

When `adata.raw` exists, expression defaults to `.raw`, matching scanpy's
plotting convention. Set `use_raw=False` for `.X` or pass `layer=` for a named
layer. A source attached to `gene(...)` overrides the plot-wide choice for that
gene.

## Main interfaces

| Interface | Purpose | Return type |
|---|---|---|
| `gganndata` | Resolve AnnData fields for the plotnine grammar | `plotnine.ggplot` |
| `plot_embedding`, `plot_features` | Embedding and feature panels | `plotnine.ggplot` |
| `plot_dotplot`, `plot_matrixplot`, `plot_heatmap` | Expression summaries and heatmaps | `plotnine.ggplot` |
| `plot_violin`, `plot_box`, `plot_ridge`, `plot_sina` | Expression distributions | `plotnine.ggplot` |
| `plot_rank_genes_dotplot`, `plot_rank_genes_matrixplot`, `plot_volcano`, `plot_ma` | Differential-expression results | `plotnine.ggplot` |
| `plot_proportions`, `plot_correlation` | Composition and correlation summaries | `plotnine.ggplot` |
| `plot_clustermap` | Clustered heatmap | `PyComplexHeatmap.ClusterMapPlotter` |
| `plot_upset` | Set intersections | marsilea UpSet object |

See the [API reference](https://mdmanurung.github.io/ggann/api.html) for all
arguments and return types. Complete helper and grammar examples are in
[`examples/grammar_equivalents.py`](examples/grammar_equivalents.py).

## Installation

```bash
pip install git+https://github.com/mdmanurung/ggann
```

Optional extras are imported only by the functions that need them:

| Extra | Functions | Backend |
|---|---|---|
| `density` | `plot_density` | [pyNebulosa](https://github.com/mdmanurung/pyNebulosa) |
| `upset` | `plot_upset` | [marsilea](https://marsilea.readthedocs.io/) |
| `heatmap` | `plot_clustermap` | [PyComplexHeatmap](https://github.com/DingWB/PyComplexHeatmap) |
| `pseudobulk` | `pseudobulk` | [decoupler](https://decoupler-py.readthedocs.io/) |

```bash
pip install "ggann[density,upset,heatmap,pseudobulk] @ git+https://github.com/mdmanurung/ggann"
```

`annplyr` and
[`plotnine-extra`](https://github.com/mdmanurung/plotnine-extra) are core
dependencies installed with ggann.

## Development

```bash
pip install -e ".[test,docs]"
pytest -q
GGANN_DOCS_OFFLINE=1 sphinx-build -W --keep-going -b html docs docs/_build/html
```

The docs build does not execute the example scripts. To regenerate committed
figures, run the relevant script under [`examples/`](examples/).
