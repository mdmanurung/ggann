# Quickstart

This guide assumes an `AnnData` named `adata` with `X_umap` in `adata.obsm`, a
`cell_type` column in `adata.obs`, and the genes used below. The
{doc}`vignettes/grammar_of_graphics` version runs during every documentation
build against a deterministic local fixture.

## One-call helper

Use a helper for a standard single-cell figure:

```python
import ggann as ag

embedding = ag.plot_embedding(
    adata,
    basis="umap",
    color="cell_type",
    label=True,
)
```

Other common one-call summaries follow the same naming conventions:

```python
genes = ["CD3D", "NKG7", "MS4A1", "CST3"]

dotplot = ag.plot_dotplot(adata, genes, group_by="cell_type")
matrix = ag.plot_matrixplot(adata, genes, group_by="cell_type")
violin = ag.plot_violin(adata, ["CD3D"], group_by="cell_type")
```

These return ordinary `plotnine.ggplot` objects.

## Grammar-style plot

Use `gganndata` when plotnine already expresses the figure clearly. Include a
field in the aesthetic mapping when a later facet or layer needs that column;
`group=` is a convenient non-visual mapping for `condition` here.

```python
from ggann import aes, gganndata, obs, obsm
from plotnine import facet_wrap, geom_point, scale_color_brewer, theme_classic

plot = (
    gganndata(
        adata,
        aes(
            x=obsm("umap", 0),
            y=obsm("umap", 1),
            color=obs("cell_type"),
            group=obs("condition"),
        ),
        add_theme=False,
    )
    + geom_point(size=1.8, alpha=0.85)
    + scale_color_brewer(type="qual", palette="Set2")
    + facet_wrap("condition")
    + theme_classic()
)
```

The returned object is a real `plotnine.ggplot`: add any compatible geom,
stat, scale, coordinate system, facet, label, or theme with `+`.

## Choose the expression source

The same source arguments apply to the grammar and expression helpers.

```python
# adata.X
x_plot = ag.plot_embedding(adata, "umap", color="CD3D", use_raw=False)

# adata.raw.X
raw_plot = ag.plot_embedding(adata, "umap", color="CD3D", use_raw=True)

# adata.layers["counts"]
counts_plot = ag.plot_embedding(adata, "umap", color="CD3D", layer="counts")
```

With neither argument, expression uses `adata.raw` when present and otherwise
`.X`, matching Scanpy's plotting convention. `layer=` and `use_raw=True` are
mutually exclusive.

Pin individual genes when one grammar plot mixes sources:

```python
from ggann import gene

mixed = gganndata(
    adata,
    aes(
        x=gene("CD3D", use_raw=True),
        y=gene("NKG7", layer="counts"),
        color=obs("cell_type"),
    ),
) + geom_point()
```

## Bound grammar materialization

`gganndata` can reject a mapping before its first matrix read when the complete
request exceeds a known boundary. The example below resolves two `obsm`
coordinates and one gene, or `3 × n_obs` logical matrix values; observation
metadata is not charged.

```python
bounded = gganndata(
    adata,
    aes(
        x=obsm("umap", 0),
        y=obsm("umap", 1),
        color=gene("CD3D", use_raw=False),
    ),
    max_matrix_values=3 * adata.n_obs,
) + geom_point()
```

An invalid or exceeded budget raises `annplyr.AnnplyrError`. High-level plotting
helpers do not currently expose this argument; use the grammar or extract a
custom bounded table with annplyr when a hard limit is required.

## Refine and save

Helper output remains composable:

```python
from plotnine import labs, theme

final = dotplot + labs(title="Lineage markers") + theme(figure_size=(7, 4))
final.save("markers.pdf", width=180, height=100, units="mm")
```

Continue with {doc}`concepts` for name resolution, ordering, missing values,
downsampling, ownership, and return types.
