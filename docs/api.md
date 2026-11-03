# API reference

```{eval-rst}
.. currentmodule:: ggann
```

## Conventions

| Concept | Rule |
|---|---|
| `genes` and `features` | `genes` selects expression variables. `features` may also select numeric observation columns where documented. Repeated sequence entries are removed in first-seen order; grouped marker maps reject a gene assigned to multiple labels. |
| `group_by` and `split_by` | Every helper that takes `group_by=` also accepts Scanpy's `groupby=`; giving both with different values raises `TypeError`. `group_by` defines the primary summary or colour groups. `split_by` adds a within-group split or facet. Cells with a missing `group_by` or `split_by` value are dropped from every grouped plot and statistic, with a warning that gives the count; proportions use the kept cells. |
| Expression source | `layer=` and `use_raw=True` are mutually exclusive. With neither set, expression uses `adata.raw` when present and otherwise `adata.X`. A source set on `gene(...)` overrides the plot-wide source. |
| Ordering | Categorical observation order is preserved. An explicit `categories_order` must include every observed, non-missing group; unused requested levels are dropped. |
| Downsampling | `downsample=` caps total cells for embedding/feature panels and cells per group for grouped distributions and heatmaps. `random_state=0` is deterministic; `None` draws a fresh sample. `random_state` also seeds the point jitter of `plot_box`, `plot_violin` and the spread of `plot_sina`, so a seeded figure has the same point positions on every draw and export; `None` is nondeterministic. |
| Materialization budget | `gganndata(..., max_matrix_values=...)` bounds cumulative expression and `obsm` values before extraction. Observation metadata is free. High-level plotting helpers do not currently expose the hard limit. |
| Colours | Categorical scales reuse `adata.uns["<column>_colors"]` when valid; otherwise they use Scanpy's default palette (20, 28 or 102 colours), and plotnine hues beyond 102 categories. `gganndata` adds no default colour scale: add `+ scale_color_obs(adata, key)` to use stored colours. For continuous colour, an explicit `cmap` wins over the style's sequential map, which wins over the family default (`cmap=None`). Numeric values use a continuous colour map. `color` is canonical; `colour` scale aliases remain available. |
| Publication mode | `style_context(...)` changes family styling and rendering defaults without changing prepared tables. Later plotnine additions win. |
| Rasterization | Dense embedding points and matrix tiles remain vectors unless `rasterized=True`; text, axes, guides, tags, and annotations remain vectors. |
| Missing colour values | In `plot_embedding`, cells with a missing `color=` value stay in the plot, drawn in `na_color` and labelled `NA` in the legend. |
| Return boundary | Grammar and plotnine-native helpers return composable `plotnine.ggplot` objects. `plot_clustermap`, `tidy_heatmap` and `plot_upset` return their grid-backend objects. |

## Performance and ownership by family

| Family | Important behavior |
|---|---|
| `gganndata`, embeddings, QC scatter | Resolve only mapped fields; per-cell rendering still scales with retained observations. Grammar calls can enforce `max_matrix_values`. |
| Dotplot, matrixplot, correlation | Project requested genes before annplyr aggregation; no implicit downsampling. |
| Violin, box, sina, ridge, tracks, heatmap | Prepare long per-cell tables. Use explicit `downsample=` where supported and scientifically appropriate. |
| `plot_highest_expr_genes` | Reads the selected whole expression source by definition, remains sparse through ranking, then materializes only top genes. |
| Density and clustering | KDE and hierarchical clustering add compute beyond extraction; optional backends document their own object and memory boundaries. |

Plot construction leaves the input `AnnData` unchanged. Prepared pandas data
are owned by the returned object; `set_theme` is the documented process-global
state change. See {doc}`grammar` and {doc}`notes` for the contracts and timings.

## Family notes

Arguments and limits that differ from what the signature alone suggests.

**Grammar.** `gganndata(adata, mapping, *, layer, use_raw, max_matrix_values, add_theme=True)`
returns a `GGAnnPlot`, a `plotnine.ggplot` subclass. Fields added later with
`+ aes(...)`, `+ geom_*(aes(...))` or `+ facet_*(...)` are read from the AnnData
when added and count against `max_matrix_values`. They do not work on the
output of `plot_*` helpers; use `split_by=` there. `add_theme=False` keeps
plotnine's default `theme_gray`. See {doc}`grammar`.

**Embeddings.**

- `plot_embedding` takes `color=[...]` (one panel per colour, `ncol=`),
  `components=`, `palette=`, `color_map=`, `vmin`/`vmax`/`vcenter` (numbers or
  `"p99"`), `groups=`, `sort_order=`, `na_color=`, `title=` and
  `legend_loc=` (`"right margin"`, `"right"`, `"on data"`, `"none"`). The
  Matplotlib backend rejects these with `ValueError`.
- `plot_embedding(method="binned")` is an opt-in approximate point density
  (2-D histogram plus Gaussian filter); the legend reads `density (binned)`.
  The exact KDE is the default and warns above 20,000 cells.
- `plot_features` draws one panel and colour bar per feature when a feature is
  an `obs` column (`color_scale="auto"`); genes share one scale
  (`color_scale="shared"` or `"free"` forces either). It rejects categorical
  features by name.
- `plot_embedding_density` adds `groups=`, `background=`, `components=`,
  `vmin=` and `vmax=`. It raises `ValueError` when `group_by` has more than
  100 levels; a numeric column counts one level per value.
- Not available in `plot_embedding`: `projection`, `edges`, `arrows`,
  `add_outline`, `frameon`, `legend_fontsize`, `legend_fontweight`, and `size`
  as an array. Not available in dotplot and matrixplot: `dendrogram`,
  `var_group_rotation`, `smallest_dot`, `mean_only_expressed`. Not available in
  `plot_embedding_density`: `key`.

**Markers and summaries.**

- The first category is drawn at the top in `plot_dotplot` and
  `plot_matrixplot`, as in Scanpy. Use `swap_axes=True`, not `coord_flip()`.
- `plot_dotplot`, `plot_matrixplot` and `plot_heatmap` take `swap_axes`,
  `vmin`/`vmax` and `cmap=None`; `plot_dotplot` also takes `dot_min`/`dot_max`.
  They warn above 100 `group_by` levels (dotplot, matrixplot).
- `genes=` may be a dict `{marker group: [genes]}`, and
  `plot_dotplot_grouped` / `plot_matrixplot_grouped` take the same mapping.
  Each marker group is one facet panel, side by side, not a Scanpy-style
  bracket above shared columns. A gene listed under several groups appears once
  in each panel, with values from one computation per gene. Plotnine backend only.
- `plot_tracksplot` accepts a `{label: [genes]}` dict; tracks are titled
  `gene [label]`.
- `plot_stacked_violin` fills each violin by group median expression
  (`cmap="Blues"`), puts groups in rows and genes along x (`swap_axes=True` transposes it), and has no y tick labels unless
  `yticklabels=True`.
- `plot_clustermap` draws no dendrogram unless you pass `row_dendrogram=True`
  or `col_dendrogram=True`. PyComplexHeatmap sizes the figure when the plot is
  constructed, so create the Matplotlib figure at the export size first (for
  example `plt.figure(figsize=(89 / 25.4, 70 / 25.4))`). Resizing afterwards
  clips dendrograms, colour bars and bottom rows. Clustering needs at least two
  entries on the axis.
- `plot_ridge` raises above `MAX_RIDGE_GROUPS = 100` groups, validates `n_grid`
  and `scale`, and takes `bw_adjust`.

**Distributions.**

- `plot_violin`, `plot_box`, `plot_sina`, `plot_stacked_violin`, `plot_ridge`
  and `plot_expression_bar` accept numeric `adata.obs` columns in `genes`.
  A bare name that is also a gene stays a gene.
- `log=True` (violin, box, sina, stacked violin) plots `log1p(expression)`;
  zeros stay on the plot. `order=` subsets and orders groups; `rotation=` turns
  the x labels.
- `plot_sina` inherits a plotnine limit: when every value in a panel is an
  integer, such as raw counts, `stat_sina` adds jitter to the heights, so
  points sit off the data values. Use a normalised layer for sina plots.
- `plot_expression_line(bins=N)` bins a numeric `x`. `plot_proportions`
  rejects `split_by` equal to `group_by`.

**Differential expression.**

- `plot_rank_genes_dotplot` and `plot_rank_genes_matrixplot` default to
  `n_genes=10`; the matrixplot defaults to `standard_scale=None`. They take
  `groups=`, `min_logfoldchange=` and `values_to_plot=` (`scores`,
  `logfoldchanges`, `pvals`, `pvals_adj`, `log10_pvals`, `log10_pvals_adj`).
  `values_to_plot` draws the stored statistic and needs a result stored with
  `n_genes=adata.n_vars`. A gene that is a top marker of several groups gets one
  column per group, labelled `GENE (group)`. `n_genes <= 0` raises
  `ValueError`; Scanpy's bottom-ranked selection is not supported.
- `plot_volcano(adata, group, ...)` needs one group name as a `str`.
- `plot_ma(data, mean=, lfc_col="log2FoldChange", pval="padj", lfc=1.0,
  padj=0.05)`: `lfc_col` names the fold-change column; `lfc` and `padj` are
  thresholds, as in `plot_volcano`. A gene is significant when `padj` is
  below the threshold and `|LFC| >= lfc`. Pass `lfc=0` to colour on `padj` alone.
  Rows with a missing or non-positive `mean` are dropped with a warning.
- `plot_upset` sorts intersections largest first; pass `ascending=True` for
  the smallest first.

**Correlation and QC.** `plot_correlation` correlates group-mean gene profiles
by default; Scanpy correlates PCA-space means, which `use_rep="X_pca"` (with
`n_pcs=`) reproduces. `plot_variance_ratio` defaults to `n_pcs=30, log=False`.
`plot_highest_expr_genes` shows 30 genes, takes `n_top` as an alias of `n`, and
adds `gene_symbols=` and `log=`. `plot_qc_violin` and `plot_qc_scatter` take
`log=`.

**Density.** `plot_density` takes `downsample=` and `random_state=`.

## Grammar

```{eval-rst}
.. autosummary::
   :toctree: generated/native
   :nosignatures:

   gganndata
   gene
   obs
   obsm
   embedding_coords
```

`ggann.aes` is the plotnine
[`aes`](https://plotnine.org/reference/aes.html) class re-exported for
convenience.

## Embeddings

```{eval-rst}
.. autosummary::
   :toctree: generated/native
   :nosignatures:

   plot_embedding
   plot_features
   plot_density
   plot_embedding_density
```

## Markers and expression summaries

```{eval-rst}
.. autosummary::
   :toctree: generated/native
   :nosignatures:

   plot_dotplot
   plot_dotplot_grouped
   plot_matrixplot
   plot_matrixplot_grouped
   plot_heatmap
   plot_violin
   plot_ridge
   plot_stacked_violin
   plot_tracksplot
   plot_clustermap
   tidy_heatmap
   TidyHeatmap
   plot_dendrogram
```

`tidy_heatmap` and its `TidyHeatmap` verbs are described in {doc}`heatmaps`.

## Distributions

```{eval-rst}
.. autosummary::
   :toctree: generated/native
   :nosignatures:

   plot_box
   plot_sina
   plot_expression_bar
   plot_expression_line
```

## Differential expression

```{eval-rst}
.. autosummary::
   :toctree: generated/native
   :nosignatures:

   rank_genes_df
   plot_rank_genes_dotplot
   plot_rank_genes_matrixplot
   plot_volcano
   plot_ma
```

## Composition, correlation, and sets

```{eval-rst}
.. autosummary::
   :toctree: generated/native
   :nosignatures:

   plot_proportions
   plot_correlation
   plot_upset
```

## Pseudobulk

```{eval-rst}
.. autosummary::
   :toctree: generated/native
   :nosignatures:

   pseudobulk
```

## Quality control

```{eval-rst}
.. autosummary::
   :toctree: generated/native
   :nosignatures:

   plot_qc_violin
   plot_qc_scatter
   plot_highest_expr_genes
   plot_variance_ratio
```

## Publication design and exact export

```{eval-rst}
.. autosummary::
   :toctree: generated/native
   :nosignatures:

   PublicationStyle
   publication_style
   theme_publication
   style_context
   publication_palette
   save_publication
```

The `single-column` and `double-column` presets are generic 89 mm and 183 mm
starting widths. They coordinate final-size typography and rendering but do not
claim compliance with a named journal. See {doc}`publication` for palette
validation, font fallback, rasterization, composition, and vector-editing
contracts.

## Exploratory scales and theme

```{eval-rst}
.. autosummary::
   :toctree: generated/native
   :nosignatures:

   theme_ggann
   set_theme
   reset_theme
   sizes
   scale_color_obs
   scale_fill_obs
   obs_colors
   scale_color_expression
   scale_fill_expression
   scale_color_celltype
   scale_fill_celltype
```

The `scale_colour_*` names are aliases of the corresponding `scale_color_*`
functions.

`set_theme(base_size=...)` is ggann's counterpart of Scanpy's
`sc.set_figure_params(fontsize=...)`; ggann has no `set_figure_params`. It is
global and applies to plots built afterwards, not to plots already built.
`style_context` overrides it inside its block and restores the previous theme,
`sizes` and rcParams on exit, so a `set_theme` call made inside the block is
discarded. `reset_theme()` installs `theme_gray` and resets `sizes` to a base
size of 11; it does not restore an earlier custom theme (use `style_context`
for that). `base_size` must be a finite positive number.

`scale_color_obs` and `scale_fill_obs` give numeric `obs` columns a continuous
scale and raise `KeyError` for a missing column. `scale_color_expression` and
`scale_fill_expression` take `cmap=None`, which means Reds (colour) or viridis
(fill) outside a style context and the style's sequential map inside one.

## Layout

```{eval-rst}
.. autosummary::
   :toctree: generated/native
   :nosignatures:

   compose
   tag_panels
```

The following names are re-exported from plotnine or
[`plotnine-extra`](https://github.com/mdmanurung/plotnine-extra). Their upstream
documentation defines their arguments and return types.

| Names | Purpose |
|---|---|
| `Wrap`, `Stack`, `Beside` | Plot composition operators |
| `plot_layout`, `plot_annotation` | Composition configuration |
| `geom_text_repel`, `geom_label_repel` | Repelled text labels |
| `ggsave` | Save a plotnine figure |
| `MatplotlibGGPlot` | Class returned by `backend="matplotlib"`; a `ggplot` subclass |

## Statistical layers

These re-exports come from
[`plotnine-extra`](https://github.com/mdmanurung/plotnine-extra):

| Names | Purpose |
|---|---|
| `stat_compare_means`, `stat_pwc`, `stat_pvalue_manual` | Group comparisons and p-value annotations |
| `stat_cor`, `stat_regline_equation` | Correlation and regression annotations |
| `stat_anova_test`, `stat_kruskal_test` | Omnibus tests |
| `stat_central_tendency` | Mean or median annotation |
| `geom_signif` | Significance brackets |

What each statistic tests and which p-value it shows is in {doc}`grammar`.
