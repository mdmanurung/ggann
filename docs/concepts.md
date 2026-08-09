# API concepts

## Two interfaces

`ggann` has one data contract and two plotting interfaces:

| Interface | Use it when | Return |
|---|---|---|
| `gganndata(adata, aes(...))` | Per-cell fields map directly to plotnine aesthetics | `plotnine.ggplot` |
| `plot_*` helpers | The plot needs aggregation, ordering, or specialized preparation | Usually `plotnine.ggplot` |

`plot_clustermap` and `plot_upset` are deliberate exceptions: they use
PyComplexHeatmap and marsilea, respectively, and return those backend objects.

## Aesthetic resolution

`gganndata` eagerly creates the minimal per-cell table required by its aesthetic
mapping. A bare string is resolved in this order:

1. a column in `adata.obs`;
2. a gene in the selected expression source;
3. a coordinate exposed by `adata.obsm`, such as `"UMAP_1"`.

An observation column wins when it shares a name with a gene. `ggann` warns
about that collision. Use a typed selector to make intent explicit:

```python
from ggann import aes, gene, gganndata, obs, obsm
from plotnine import geom_point

plot = gganndata(
    adata,
    aes(
        x=obsm("umap", 0),
        y=obsm("umap", 1),
        color=gene("CD3D", use_raw=True),
        shape=obs("condition"),
    ),
) + geom_point()
```

The equivalent string forms are useful in configuration files:

| Selector | Source |
|---|---|
| `"obs:condition"` | `adata.obs["condition"]` |
| `"gene:CD3D"` | the plot-wide expression source |
| `"gene:CD3D@X"` | `adata.X` |
| `"gene:CD3D@raw"` | `adata.raw.X` |
| `"gene:CD3D@counts"` | `adata.layers["counts"]` |
| `"obsm:umap[0]"` | zero-based coordinate 0 of `X_umap` |

Typed selectors raise source-specific errors when a requested source is absent.
They are preferable to bare strings in reusable code.

Only fields referenced by `aes(...)` are added to the plot data. If a later
facet or layer needs another observation column, include it in the mapping—for
example `group=obs("condition")` before adding `facet_wrap("condition")`.

## `.X`, `.raw`, and layers

Expression-source selection is consistent across grammar and helpers:

| Arguments | Source |
|---|---|
| `use_raw=False` | `adata.X` |
| `use_raw=True` | `adata.raw.X` |
| `layer="counts"` | `adata.layers["counts"]` |
| neither | `.raw` when present, otherwise `.X` |

`layer=` and `use_raw=True` cannot be combined. A `gene(...)` selector with its
own source overrides the plot-wide choice for that gene, which permits
intentional mixed-source grammar plots.

## Grouping and splitting

`group_by` defines the primary discrete group used for summaries, colours, or
the distribution x-axis. `split_by` introduces a second within-group split or a
facet where the helper supports it.

- Categorical `obs` order is preserved.
- Non-categorical groups use sorted observed values when they are mutually
  comparable; otherwise first-seen order is retained.
- `categories_order=` overrides the primary order. It must contain every
  observed non-missing group; unused requested levels are dropped.
- Grouped summaries omit missing grouping values. Per-cell plots keep them in
  the prepared table and plotnine applies the selected geom's missing-value
  behavior.

Use a pandas categorical when one order should apply across every figure:

```python
adata.obs["cell_type"] = adata.obs["cell_type"].cat.reorder_categories(
    ["T cell", "NK cell", "B cell", "Monocyte"]
)
```

Use `categories_order=` for a plot-local order:

```python
ag.plot_dotplot(
    adata,
    ["CD3D", "NKG7", "MS4A1", "CST3"],
    group_by="cell_type",
    categories_order=["T cell", "NK cell", "B cell", "Monocyte"],
)
```

## Downsampling

`downsample=` is an explicit opt-in, never an implicit performance shortcut.

- Embedding and feature helpers cap the total number of retained cells.
- Grouped distribution and per-cell heatmap helpers cap cells within each
  group.
- `random_state=0` gives deterministic sampling. Pass another integer for a
  different reproducible sample, or `None` for a fresh sample.
- A non-positive or non-integer `downsample` value raises `ValueError`.

Aggregated helpers such as dotplot and matrixplot summarize all cells unless
their documented API says otherwise.

## annplyr projection and materialization

Tabular extraction is delegated to annplyr. Before conversion, ggann projects
expression to the requested genes and observation columns; sparse inputs stay
sparse through that bounded projection and are densified only when plotnine
needs its final small table.

`gganndata(..., max_matrix_values=...)` preflights cumulative expression and
`obsm` values across the complete mapping before its first accessor read.
Observation metadata is not charged. `embedding_coords` exposes the same
boundary for direct coordinate extraction. A negative or exceeded budget raises
the public `annplyr.AnnplyrError`.

```python
bounded_plot = gganndata(
    adata,
    aes(
        x=obsm("umap", 0),
        y=obsm("umap", 1),
        color=gene("CD3D", use_raw=False),
    ),
    max_matrix_values=3 * adata.n_obs,
) + geom_point()
```

High-level plotting helpers bound *width* by requested features but do not
currently expose the hard budget. Use the grammar, or extract explicitly with
annplyr and plot the resulting table:

```python
table = adata.ap.to_df(
    obs=["cell_type"],
    x=["CD3D", "NKG7"],
    max_matrix_values=2 * adata.n_obs,
)
```

The annplyr planner validates all requested matrix sources before its first
read. See the executable {doc}`vignettes/annplyr_interop` workflow.

Dense, CSR, CSC, views, duplicate observation names, and backed inputs are
supported by the projected paths covered in the compatibility suite. Backed
plotting necessarily materializes the requested plotting table in memory; it
does not load unrelated expression columns.

## Ownership and mutation

Plot construction does not mutate the input `AnnData`. Downsampling uses a
view, expression extraction delegates projected public annplyr frames, and
prepared pandas tables are owned by the returned plot. Duplicate observation
names are handled positionally during projection and restored in prepared
output.

Two operations affect state outside `AnnData`:

- `set_theme(...)` changes plotnine's process-wide default theme;
- saving or drawing a plot creates renderer state and output files.

Call `reset_theme()` to restore plotnine's default theme.

## Return and composition boundary

| Functions | Return |
|---|---|
| `gganndata`, plotnine-native `plot_*` helpers | `plotnine.ggplot` |
| `compose` | plotnine composition object |
| `embedding_coords`, `rank_genes_df` | `pandas.DataFrame` |
| `plot_clustermap` | `PyComplexHeatmap.ClusterMapPlotter` |
| `plot_upset` | marsilea UpSet object |
| `pseudobulk` | `AnnData` |

Ordinary plots compose with plotnine's `+`. Combine panels with `compose(...)`
or plotnine's composition operators, then save at an exact physical size.
