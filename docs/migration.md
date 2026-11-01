# Coming from Scanpy or scplotter

`ggann` is useful when an AnnData plotting workflow needs reusable plot objects,
plotnine composition, or one coherent path from exploration to exact-size
publication export. The data model remains AnnData; there is no conversion to a
second analysis object.

## The main mental-model shift

Scanpy plotting functions usually draw immediately. `ggann` prepares the same
kind of biological summary and returns a plot object:

```python
# Scanpy
sc.pl.dotplot(adata, genes, groupby="louvain", show=False)

# ggann
dotplot = ag.plot_dotplot(adata, genes, group_by="louvain")
```

You can display `dotplot`, add plotnine components with `+`, place it in a
composition, or export it later. Except for clustermap and UpSet, helper output
is an ordinary `plotnine.ggplot`.

## Naming map

| Task | ggann spelling | Practical meaning |
|---|---|---|
| primary category | `group_by` (`groupby=` also accepted) | category used for summaries, fills, or rows |
| secondary category | `split_by` | facet or within-group split |
| expression source | `layer` / `use_raw` | named layer, `.raw`, or `.X` |
| category order | `categories_order` | complete order of observed non-missing groups |
| explicit cell cap | `downsample` | never applied silently |
| reproducible sampling | `random_state` | deterministic `0` by default |
| grammar read bound | `max_matrix_values` | cumulative gene and `obsm` values before extraction |
| colour aesthetic | `color` | canonical spelling; `scale_colour_*` aliases are available |

## Common Scanpy translations

The examples use the 2,638-cell PBMC3K dataset, which Scanpy downloads once
and then caches. Highly-variable-gene selection dropped `CD3D` and `LYZ` from
`.X`, so marker calls read `.raw` with `use_raw=True`:

```python
import numpy as np
import scanpy as sc
import ggann as ag

adata = sc.datasets.pbmc3k_processed()
genes = ["CD3D", "MS4A1", "NKG7", "GNLY", "CST3"]

# A median split of the measured library size, used by the facets below.
adata.obs["depth"] = np.where(
    adata.obs["n_counts"] < adata.obs["n_counts"].median(), "low", "high"
)
```

### UMAP by annotation

```python
# Scanpy
sc.pl.umap(adata, color="louvain", legend_loc="on data", show=False)

# ggann
umap = ag.plot_embedding(
    adata,
    "umap",
    color="louvain",
    label=True,
)
```

The ggann version keeps both direct labels and the legend. Use
`rasterized=True` for a dense vector figure, or `downsample=` only when the
scientific use permits an explicit subset.

### Marker dotplot and matrixplot

```python
# Scanpy
sc.pl.dotplot(adata, genes, groupby="louvain", use_raw=True, show=False)
sc.pl.matrixplot(adata, genes, groupby="louvain", use_raw=True, show=False)

# ggann
dotplot = ag.plot_dotplot(
    adata, genes, group_by="louvain", use_raw=True
)
matrix = ag.plot_matrixplot(
    adata, genes, group_by="louvain", use_raw=True
)
```

Dot colour is mean expression and size is fraction above
`expression_cutoff`. Matrix tiles are raw group means, as in Scanpy
(`standard_scale=None`); pass `standard_scale="var"` to rescale each gene to
`[0, 1]` in either library. Both plots draw the first category at the top, as
Scanpy does. To put genes on y and groups on x, use `swap_axes=True`, not
`coord_flip()`, which reverses the displayed order. In publication mode,
`annotate="auto"` adds contrast-aware values only when a final-size cell is at
least 12 points wide and high.

### Violin with a visible centre and interval

```python
# Scanpy
sc.pl.violin(adata, "CD3D", groupby="louvain", use_raw=True, show=False)

# ggann
violin = ag.plot_violin(
    adata,
    ["CD3D"],
    group_by="louvain",
    use_raw=True,
    add_box=True,     # median and interquartile range
    add_points=False,
    stats=False,
)
```

The KDE and inner box use every cell unless `downsample` is set. If
`stats=True`, the test is computed on the retained cells, so leave downsampling
off when the p-value must represent the full population.

### Cell-type composition

```python
composition = ag.plot_proportions(
    adata,
    group_by="louvain",
    split_by="depth",
    normalize=True,
)
```

Each depth bar sums to one. The helper counts observation metadata only and
does not materialize expression.

### Stacked violin

```python
# Scanpy
sc.pl.stacked_violin(
    adata, genes, groupby="louvain", use_raw=True, swap_axes=True, show=False
)

# ggann
stacked = ag.plot_stacked_violin(adata, genes, group_by="louvain", use_raw=True)
```

ggann draws genes as rows and groups along x. Scanpy 1.12 draws the transpose by
default, so pass `swap_axes=True` to Scanpy for the same layout. Each violin is
filled by the group's median expression with one shared colour bar
(`cmap="Blues"`), as in Scanpy. Rows have no y tick labels; pass
`yticklabels=True` to show them.

### Repeated markers

A gene that belongs to several marker groups is drawn once per group in
`plot_dotplot`, `plot_matrixplot`, `plot_heatmap` and `plot_tracksplot`:

```python
marker_sets = {"T": ["CD3D", "NKG7"], "NK": ["NKG7", "GNLY"]}
repeated = ag.plot_dotplot(adata, marker_sets, group_by="louvain", use_raw=True)
```

The panels are labelled by marker group, and `NKG7` appears in both with values
from one computation. Scanpy draws brackets instead; ggann does not.

### Differential-expression views

```python
sc.tl.rank_genes_groups(adata, "louvain", method="wilcoxon", key_added="de")

# Scanpy
sc.pl.rank_genes_groups_dotplot(adata, n_genes=5, key="de", show=False)

# ggann
ranked = ag.plot_rank_genes_dotplot(adata, n_genes=5, key="de", groupby="louvain")
```

`n_genes` defaults to 10, and `plot_rank_genes_matrixplot` uses raw group means
(`standard_scale=None`), as Scanpy does. Negative `n_genes` is not supported.

For a pseudobulk MA plot, `lfc_col` names the fold-change column and `lfc` is a
threshold, as in `plot_volcano`:

```python
ag.plot_ma(table, mean="baseMean", lfc_col="log2FoldChange", lfc=1.0, padj=0.05)
```

Earlier builds read `lfc=` as the column name. Pass `lfc=0` to colour on `padj`
alone.

## scplotter-style ergonomic options

The helper vocabulary includes direct embedding labels, `split_by` facets,
grouped marker maps, nested violin boxes, sina/ridge distributions, expression
summaries, differential-expression views, composition plots, correlations,
dendrograms, clustermaps, and UpSet diagrams. These stay in Python/AnnData and
return plotnine objects wherever the underlying geometry permits it.

For a custom layout, compose helpers rather than configuring a bespoke plotting
class:

```python
figure = ag.compose(
    [umap, dotplot, violin, composition],
    ncol=2,
    widths=(0.95, 1.25),
    heights=(1.0, 1.0),
    gap=2,
    tag_levels="auto",
)
```

## When to use the grammar directly

Use `gganndata` when plotnine already describes the figure:

```python
from ggann import aes, gene, gganndata, obs, obsm
from plotnine import facet_wrap, geom_point

depth_umap = (
    gganndata(
        adata,
        aes(
            x=obsm("umap", 0),
            y=obsm("umap", 1),
            color=gene("CD3D", use_raw=True),
            group=obs("depth"),
        ),
    )
    + geom_point(size=1)
    + facet_wrap("depth")
)
```

Explicit selectors remove ambiguity when an observation column and a gene have
the same name. Bare names resolve through observation metadata, the selected
expression source, then embedding coordinates.

## Expression-source rules

With neither `layer` nor `use_raw` set, ggann uses `.raw` when available and
otherwise `.X`, matching Scanpy's plotting convention. Use `use_raw=False` to
force `.X`; a named `layer` selects that layer. Combining `layer=` with
`use_raw=True` is an error. A source attached to `gene(...)` overrides the
plot-wide selection for that gene.

Only requested variables are projected through annplyr before conversion or
aggregation. Grammar calls can set `max_matrix_values` to reject an oversized
gene/embedding request before the first matrix read.

## Category, palette, and missing-value behavior

Categorical order in `adata.obs` is preserved. An explicit
`categories_order` must cover every observed, non-missing group. Grouped
summaries omit rows with missing grouping keys.

Categorical plots reuse `adata.uns["<obs>_colors"]`. Publication mode validates
exact category alignment and colour syntax; an invalid palette warns and falls
back deterministically without changing the AnnData. Missing observations in
embeddings remain present and use neutral grey.

## Move the assembled figure to publication mode

The calls do not need to be rewritten:

```python
with ag.style_context("double-column"):
    panels = [
        ag.plot_embedding(adata, "umap", color="louvain", label=True),
        ag.plot_dotplot(adata, genes, group_by="louvain", use_raw=True),
    ]
    figure = ag.compose(panels, ncol=2, widths=(0.9, 1.3), gap=2)

ag.save_publication(
    figure,
    "pbmc_figure",
    width="double-column",
    height=90,
    formats=("svg", "pdf", "png"),
    dpi=600,
)
```

Scanpy's `show=`, `return_fig=`, `save=` and `ax=` have no ggann counterpart.
Display the returned plot, call `plot.draw()`, use `ag.save_publication`, or
place plots in a grid with `ag.compose`. {doc}`scanpy_parity` lists every
supported call and keyword.

See {doc}`vignettes/scanpy_migration` for an executable translation,
{doc}`scanpy_parity` and {doc}`scplotter_parity` for detailed coverage, and
{doc}`publication` for final-size design and export contracts.
