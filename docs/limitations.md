# Known limitations

This page lists what does not work in ggann 0.1.0 and what to do instead. Each
line states the limit, then the workaround. Fixed items are not listed. Items
marked "upstream" come from plotnine or plotnine-extra and are reproduced in
plain plotnine code; none has been filed yet.

## Embeddings

- **3-D embeddings.** `plot_embedding` has no `projection="3d"`. Plot two
  components with `components=`, or use Matplotlib directly.
- **Missing scanpy options.** `plot_embedding` has no `edges`, `arrows`,
  `add_outline`, `frameon`, `legend_fontsize`, `legend_fontweight`, or an array
  for `size`. Add plotnine layers or a `theme()` instead.
- **Matplotlib backend options.** `backend="matplotlib"` raises `ValueError` for
  `components`, `palette`, `color_map`, `vmin`, `vmax`, `vcenter`, `groups`,
  `sort_order`, `na_color`, `title` and `legend_loc`, for dict `genes`, and for
  `split_by`. Use the plotnine backend for these.
- **Facets on a Matplotlib plot.** `plot_embedding(backend="matplotlib") +
  facet_wrap(...)` fails because the plot data holds only x, y and the colour.
  Use `split_by=` with the plotnine backend.
- **Exact point density is slow.** The default KDE took 27 s at 20,000 cells and
  warns above that size. Use `method="binned"` (approximate, labelled in the
  legend), `pointdensity=False` or `downsample=`. The speed target set for the
  exact method is unmet.
- **`plot_embedding_density` has no `key=`.** Densities are computed in ggann, so
  a stored `scanpy.tl.embedding_density` result cannot be passed in.
- **`plot_features` return type.** With an obs column among the features it
  returns a composition, not one `ggplot`. Pass `color_scale="shared"` to keep
  one plot, but a shared scale mixes units.

## Dotplot, matrixplot, heatmap and tracksplot

- **No `dendrogram=`.** `plot_dotplot`, `plot_matrixplot` and `plot_heatmap` do
  not order or draw a dendrogram. Set the group order with `categories_order=`, or use
  `plot_clustermap` or `plot_dendrogram`.
- **Missing scanpy options.** `var_group_rotation`, `smallest_dot` and
  `mean_only_expressed` do not exist. Use `dot_min` and `dot_max` for the dot
  size range.
- **No `expression_cutoff` on matrixplot and heatmap.** They show means, so
  there is no detection fraction to cut. Use `plot_dotplot`.
- **Dict `genes` needs the plotnine backend.** The Matplotlib backend raises for
  `{marker group: [genes]}`.
- **`coord_flip()` reverses the category order.** Use `swap_axes=True`.

## Rank-genes and differential expression

- **Negative `n_genes`.** Scanpy's bottom-ranked selection is not supported.
  `n_genes <= 0`, non-integers and `None` raise `ValueError`. Sort the DE table
  yourself and pass the genes to `plot_dotplot`.
- **No marker brackets.** A gene that marks several groups gets one labelled
  column per group, not a bracket over a block.
- **`values_to_plot` needs a full DE result.** Run `rank_genes_groups` with
  `n_genes=adata.n_vars`; the statistic is never an expression mean and cannot
  be combined with `standard_scale`.
- **`plot_correlation` statistic.** The default correlates group-mean gene
  profiles. Scanpy correlates PCA-space means; pass `use_rep="X_pca"` for that.
- **`plot_qc_violin` has no `use_raw`.** QC metrics are `obs` columns, which
  never live in `.raw`.
- **`plot_dendrogram` labels.** Leaf label rotation and colours are not
  configurable. Add a `theme()`.
- **UpSet axis titles.** `plot_upset` adds no set-size title; Marsilea draws the
  counts. Annotate in the caption.
- **Wilcoxon p-values.** The `stat_compare_means` label is SciPy's asymptotic p
  with continuity correction. R's exact test (n < 50) differs in the third digit
  (0.8620 against 0.8604 on the test fixture). State the test in the legend.

## Distributions

- **Sina on integer data (upstream, U5).** plotnine's `stat_sina` jitters y when
  every value in a panel is an integer, so points draw at heights not in the
  data (0 to 10 plotted as -0.17 to 10.09). ggann passes the data unchanged.
  Use `plot_violin` or `plot_box` for counts, or a continuous (for example
  log-normalized) value.
- **`stat_compare_means(comparisons=)` (upstream, U4).** It fails with
  `KeyError: 'PANEL'` in plain plotnine. Use the global test, or `ggann.stat_pwc`
  for pairwise labels.
- **`geom_signif(comparisons=)` (upstream, U1).** Not accepted. Use numeric
  `xmin`/`xmax` with `y_position`, or a data-frame mapping. ggann does not call
  `geom_signif`.
- **Ridge limits.** `plot_ridge` raises above 100 groups (`MAX_RIDGE_GROUPS`).
  The trim control is not offered; set the bandwidth with `bw_adjust`.
- **QC tick labels.** `plot_qc_violin` hides x tick labels with `group_by`. Its
  docstring shows how to restore them.
- **Repel labels (upstream, U6).** `geom_text_repel` takes `max_iter` only, not
  ggrepel's `max_time` or `max_overlaps`. Bound the work with `max_iter=`.

## Grammar of graphics

- **No late fields on helper outputs.** `gganndata` plots resolve fields named
  after creation (`+ aes(...)`, `+ facet_*(...)`). The output of
  `plot_embedding`, `plot_features` and other helpers does not, because its data
  is sorted, long or aggregated. Use `split_by=` or the helper's own options.
- **Late fields are a snapshot.** Fields are read when added, share
  `max_matrix_values`, and raise `ValueError` if cells were removed or reordered
  in between. Edits to values are not detected.
- **G4 leftovers.** `gganndata` has no policy for unused categories, does not
  validate a wrong `adata` or mapping type, and has no policy for duplicate
  variable names. Drop unused categories and make names unique before plotting.
- **No default colour scale in the raw grammar.** Add
  `+ scale_color_obs(adata, key)` to use stored `adata.uns[f"{key}_colors"]`.
  A default scale would print plotnine's "Scale is already present" message on
  every override.
- **Composition save (upstream, U2, U3).** `plotnine.ggsave` fails on
  compositions and `plotnine_extra.Wrap.save` fails with `no attribute 'axs'`;
  `Wrap(ncol=)` is not an argument. Use `ggann.compose(...).save(...)` or
  `save_publication`.
- **Panel tags (upstream, U7).** `plot_annotation` has no `tag_levels`. Use
  `compose(tag_levels=)`.

## Publication export and style

- **Clustermap size.** A `plot_clustermap` figure must be created at the export
  size first. Even then, long column names or a legend wider than the canvas are
  cut off; `save_publication` warns and does not fix it. Shorten the labels, or
  pass `show_colnames=False` and `plot_legend=False`.
- **Clustermap and compose.** `compose` rejects a clustermap or UpSet result
  (`TypeError`). Export those on their own and combine the files outside ggann.
- **Narrow Matplotlib dotplots.** The legend fits an 89 mm canvas on
  `pbmc68k_reduced`, but its fixed-position legend and colourbar still overflowed
  a 70 mm canvas by about 2 mm. `save_publication` reports it. Widen the canvas
  or use the plotnine backend.
- **Overflow is reported, not fixed.** `save_publication` never crops or enlarges
  the canvas.
- **Eight CVD-measured colours.** Only the first eight qualitative colours are
  measured under colour-vision simulation. Larger palettes keep earlier colours
  fixed, but separation is not guaranteed; ggann warns from 10 categories. Add
  labels or shapes.
- **Fonts.** Font fallback was tested with a simulated Arial (DejaVu Sans
  registered under that name), because the test host has no Arial. Check a real
  Arial install yourself.
- **`style_context` and threads.** The active style is not inherited by new
  threads, and concurrent use is unverified. Export from one thread.
- **Cmap divergence.** An explicit `cmap` with `standard_scale="zscore"` inside
  `style_context` is zero-centred in the plotnine backend, but the Matplotlib
  backend draws it without the zero-centred norm.
- **`reset_theme` is not an undo.** It installs `theme_gray`; it does not restore
  an earlier custom theme. Use `style_context`.
- **Export formats.** Only SVG, PDF, PNG and TIFF are written; `.jpg` and
  `.eps` raise.

## Scplotter and Scanpy parity

- **Lineage, velocity, neighbour graph and embedded-chart features** of
  scplotter are outside ggann's data model; see {doc}`scplotter_parity`.
- **Correlation pair matrices** and lineage-aware feature plots are not
  implemented.
