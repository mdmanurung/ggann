# API reference

```{eval-rst}
.. currentmodule:: ggann
```

## Conventions

| Concept | Rule |
|---|---|
| `genes` and `features` | `genes` selects expression variables. `features` may also select numeric observation columns where documented. Repeated sequence entries are removed in first-seen order; grouped marker maps reject a gene assigned to multiple labels. |
| `group_by` and `split_by` | `group_by` defines the primary summary or colour groups. `split_by` adds a within-group split or facet. Missing group values are omitted from grouped summaries. |
| Expression source | `layer=` and `use_raw=True` are mutually exclusive. With neither set, expression uses `adata.raw` when present and otherwise `adata.X`. A source set on `gene(...)` overrides the plot-wide source. |
| Ordering | Categorical observation order is preserved. An explicit `categories_order` must include every observed, non-missing group; unused requested levels are dropped. |
| Downsampling | `downsample=` caps total cells for embedding/feature panels and cells per group for grouped distributions and heatmaps. `random_state=0` is deterministic; `None` draws a fresh sample. |
| Materialization budget | `gganndata(..., max_matrix_values=...)` bounds cumulative expression and `obsm` values before extraction. Observation metadata is free. High-level plotting helpers do not currently expose the hard limit. |
| Colours | Categorical scales reuse `adata.uns["<column>_colors"]` when valid and otherwise use a qualitative scale. Numeric values use a continuous colour map. `color` is canonical; `colour` scale aliases remain available. |
| Publication mode | `style_context(...)` changes family styling and rendering defaults without changing prepared tables. Later plotnine additions win. |
| Rasterization | Dense embedding points and matrix tiles remain vectors unless `rasterized=True`; text, axes, guides, tags, and annotations remain vectors. |
| Return boundary | Grammar and plotnine-native helpers return composable `plotnine.ggplot` objects. `plot_clustermap` and `plot_upset` return their grid-backend objects. |

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
state change. See {doc}`concepts` and {doc}`performance` for the full contracts.

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
   plot_dendrogram
```

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
