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

Only fields referenced by `aes(...)` are added to the plot data.

### Late fields

`gganndata` returns a `GGAnnPlot`, a `plotnine.ggplot` subclass. Fields named
after construction are read from the same AnnData, so a facet or layer no
longer needs the column in the first mapping:

```python
from plotnine import facet_wrap, geom_point

plot = (
    gganndata(adata, aes(x=obsm("umap", 0), y=obsm("umap", 1)))
    + geom_point(aes(color="louvain"))
    + facet_wrap("louvain")
)
```

- A bare `+ aes(...)`, a layer with its own `aes` and no `data=`, and a facet
  request fields. Themes, scales, labels and layers with their own `data=`
  request nothing.
- Fields are read when the component is added, not when the plot was created.
  Editing values in the AnnData in between is not detected. Removing or
  reordering cells raises `ValueError`.
- Late fields count against the original `max_matrix_values` budget; the
  budget is shared, not renewed. Observation columns are free.
- `+` returns a new plot; the original and the AnnData are unchanged.
- This works only on `gganndata` plots. The data of `plot_*` helper outputs is
  sorted, long or aggregated, so a per-cell field cannot be joined to it. Use
  the helper's `split_by=` instead.

### Expressions in `aes`

An `aes` string may be an expression. ggann parses it with `ast` and never
evaluates it; plotnine evaluates it against the plot data. The supported subset
is names, numbers, strings, operators, calls and attribute access, for example
`"np.log1p(n_counts) / 2"` or `"factor(louvain)"`.

- A name is a field unless it is called (`factor(...)`) or is the base of an
  attribute (`np.log1p`). `np` works when the calling module imported it.
- A whole string that is itself a field is kept whole: `"CD3D"`,
  `"HLA-DRA"`, `"gene:CD3D"`.
- A gene whose name is not a Python identifier is written `Q("HLA-DRA")` inside
  an expression, or `gene("HLA-DRA")` as the whole value:

```python
import numpy as np

plot = gganndata(
    adata,
    aes("louvain", 'np.log1p(Q("HLA-DRA"))'),
)
```

Any column of any `obsm` array resolves by name (`PC_5`) or with
`obsm("X_pca", 4)`.

### Colours in the grammar

`gganndata` adds no default colour scale. Add `+ scale_color_obs(adata, key)`
to use the colours stored in `adata.uns["<key>_colors"]`. A default scale would
make plotnine print "Scale for 'color' is already present" every time you add
your own. The `plot_*` helpers apply stored palettes themselves.

Without stored colours, `scale_color_obs` and the helpers use Scanpy's default
palette (20, 28 or 102 colours) and plotnine hues beyond that. `adata.uns` is
never modified.

For continuous colour the order is: an explicit `cmap`, then the style's
sequential map inside `style_context`, then the family default. The helpers take
`cmap=None` to mean "use that order". `plot_embedding`'s `low`/`high` drive only
the legacy two-colour gradient and are ignored inside a style context.

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
- Cells whose `group_by` or `split_by` value is missing are dropped before
  aggregation, plotting and statistics, in grouped summaries and per-cell
  plots alike. One `UserWarning` names the field and the dropped count;
  nothing is printed when no cell is dropped. Proportions use the kept cells
  as the denominator. An embedding's `color=` field is not a grouping: its
  missing cells stay in the plot under the `NA` label, as in Scanpy.

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
  different reproducible sample, or `None` for a fresh sample. In `plot_box`,
  `plot_violin` and `plot_sina` the same seed (0 included) also fixes the point
  jitter and spread, so every draw and export places points identically.
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

## Limits of the helpers

- `plot_ridge` raises above `MAX_RIDGE_GROUPS = 100` groups.
  `plot_embedding_density` raises, and `plot_dotplot` and `plot_matrixplot`
  warn, when `group_by` has more than 100 levels. A numeric column counts one
  level per distinct value.
- `plot_correlation` correlates group-mean gene profiles by default. Scanpy
  correlates PCA-space group means; `use_rep="X_pca"` reproduces that.
- `plot_sina` on integer-valued data (raw counts) shows heights that are not
  the data: plotnine's `stat_sina` adds jitter when every value in a panel is an
  integer. Plot a normalised layer instead.

## Ownership and mutation

Plot construction does not mutate the input `AnnData`. Downsampling uses a
view, expression extraction delegates projected public annplyr frames, and
prepared pandas tables are owned by the returned plot. Duplicate observation
names are handled positionally during projection and restored in prepared
output.

Two operations affect state outside `AnnData`:

- `set_theme(...)` changes plotnine's process-wide default theme and
  `ggann.sizes`;
- saving or drawing a plot creates renderer state and output files.

`set_theme(base_size=...)` is the counterpart of Scanpy's
`sc.set_figure_params(fontsize=...)`; ggann has no `set_figure_params`. It
applies to plots built afterwards. `style_context` overrides it inside a block
and restores the previous theme, `sizes` and rcParams on exit.

`reset_theme()` installs `theme_gray` and resets `sizes`. It does not restore an
earlier custom theme; use `style_context` for that.

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
