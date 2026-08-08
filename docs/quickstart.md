# Quickstart

## Load example data

```python
import scanpy as sc
import ggann as ag
from ggann import aes, gene, gganndata
from plotnine import geom_point

adata = sc.datasets.pbmc68k_reduced()
markers = ["CD3D", "NKG7", "CST3"]
group = "bulk_labels"
```

## Use the grammar

`gganndata` resolves observation columns, genes, and embedding coordinates into
the DataFrame used by plotnine:

```python
plot = gganndata(
    adata,
    aes("UMAP_1", "UMAP_2", color=group),
) + geom_point(size=1.5)
```

The result is a `plotnine.ggplot`. Add plotnine layers, scales, facets, and
themes with `+`.

Use an accessor when the source is ambiguous or a specific expression matrix is
required:

```python
plot = gganndata(
    adata,
    aes("UMAP_1", "UMAP_2", color=gene("CD3D", use_raw=True)),
) + geom_point(size=1.5)
```

Bare gene names use `adata.raw` when it exists. Pass `use_raw=False` for
`adata.X`, or `layer="counts"` for a named layer.

## Use a helper

Plotnine-native helpers prepare the data and return a plotnine object:

```python
embedding = ag.plot_embedding(adata, "umap", color=group, label=True)
dotplot = ag.plot_dotplot(adata, markers, group_by=group)
violin = ag.plot_violin(adata, markers[:1], group_by=group)
heatmap = ag.plot_heatmap(adata, markers, group_by=group, standard_scale="var")
```

For `plot_heatmap`, `standard_scale` may be `None`, `"var"`, `"group"`, or
`"zscore"`.

The helper result remains composable:

```python
from plotnine import labs, theme

dotplot + labs(title="Marker expression") + theme(figure_size=(7, 4))
```

Complete helper and grammar constructions are in
[`examples/grammar_equivalents.py`](https://github.com/mdmanurung/ggann/blob/main/examples/grammar_equivalents.py).

## Set a session-wide theme

`set_theme` changes plotnine's global default theme. Call `reset_theme` to undo
that change.

```python
ag.set_theme(base_size=9, family="Arial")
themed = ag.plot_embedding(adata, "umap", color=group)
ag.reset_theme()
```

`ag.sizes` exposes the font sizes used by the current ggann theme. Use
`ag.sizes.geom(...)` when passing a point size to `geom_text`.

## Compose panels

```python
figure = ag.compose(
    [
        ag.plot_embedding(adata, "umap", color=group),
        ag.plot_dotplot(adata, markers, group),
        ag.plot_violin(adata, markers[:1], group),
        ag.plot_proportions(adata, group, split_by="phase"),
    ],
    ncol=2,
)

figure.save("figure1.pdf", width=180, height=140, units="mm")
```

Set `tag_levels` to `"a"`, `"1"`, `"i"`, or `None` to change panel labels.

## Grid-based plots

`plot_clustermap` and `plot_upset` do not return plotnine objects. Install their
optional extras and use the save methods provided by their backends. See
{doc}`installation` for the corresponding extras.
