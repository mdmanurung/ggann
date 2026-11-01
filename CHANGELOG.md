# Changelog

All notable changes to ggann are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and releases use
[semantic versioning](https://semver.org/).

## [0.1.0] - Unreleased

Initial public release.

### AnnData grammar

- A plotnine-native grammar for resolving observation fields, genes, layers,
  raw expression, and embedding coordinates through annplyr.
- Explicit `gene`, `obs`, and `obsm` selectors for unambiguous mixed-source
  plots.
- Optional cumulative `max_matrix_values` budgets on `gganndata` and
  `embedding_coords`, checked before matrix reads.
- Sparse, backed, category-order, missing-value, and prepared-data ownership
  contracts.

### Single-cell plot families

- Thirty public plotting helpers covering embeddings and density; dotplots,
  matrixplots, heatmaps, and correlations; distributions and expression
  summaries; markers and differential expression; quality control;
  proportions, dendrograms, clustermaps, and UpSet diagrams.
- Plotnine return values for twenty-eight helpers, with PyComplexHeatmap and
  Marsilea objects retained for the two documented grid-backend exceptions.
- An opt-in `MatplotlibGGPlot` path for supported high-volume embedding and
  matrix workloads while preserving ordinary plotnine fallback behavior after
  grammar additions.
- Explicit point/tile rasterization and contrast-aware `"auto"` or `"force"`
  matrix annotations.

### Publication design system

- Immutable, validated `PublicationStyle` settings with generic 89 mm
  single-column and 183 mm double-column presets.
- `theme_publication` family treatments and a nestable `style_context` that
  restores plotnine theme state, ggann sizes, Matplotlib parameters, and
  enclosing context state exactly.
- Deterministic qualitative, sequential, and zero-centred diverging palettes;
  strict AnnData palette validation; neutral missing-value representation; and
  normal/grayscale/colour-vision accessibility checks.
- Weighted plotnine-native composition with physical millimetre gaps, exact
  duplicate-guide collection, and automatic lowercase publication tags.
- `save_publication` for exact-size SVG, PDF, PNG, and TIFF output, editable SVG
  text, Type 42 PDF font paths, 300/600-DPI raster output, validated backgrounds,
  and export of the two grid-backend result types.

### Reproducibility and documentation

- Executable, offline documentation vignettes spanning Scanpy translation,
  grammar composition, sparse/backed inputs, annplyr interoperability,
  publication figures, and evidence-backed Scanpy trade-offs.
- A real `pbmc68k_reduced` first-use workflow that turns a familiar Scanpy
  marker review into reusable UMAP, dotplot, violin, and unequal-layout objects.
- A real `pbmc68k_reduced` publication workflow with a scientific claim, panel
  map, one shared colour vocabulary, stated populations and summaries, unequal
  panels, exact multi-format export, AnnData fingerprinting, and a
  machine-readable manifest.
- Reproducible performance, hotspot, visual, accessibility, packaging, and
  downstream-artifact invariance tooling.
- Ruff, Pyright, API documentation, executable vignette, offline Sphinx, and
  multi-version Python release checks.

[0.1.0]: https://github.com/mdmanurung/ggann/releases/tag/v0.1.0

### Data-handling contracts

- Cells whose `group_by` or `split_by` value is missing are dropped from every
  grouped plot, statistic and `pseudobulk` call, with one `UserWarning` giving
  the dropped count. Proportions use the kept cells as the denominator.
  Previously per-cell plots kept them and heatmaps crashed at draw.
- A bare string gene, feature or QC metric is one name; it was split into
  characters. Empty requests raise `ValueError` at the call.
- Helper-requested fields that resolve to nothing raise `KeyError` at the
  call instead of being skipped or failing at draw; the error names the other
  expression source when the gene lives there.
- `plot_features` rejects categorical features by name instead of dropping
  them silently.
- Boolean observation columns colour embeddings with a discrete scale.

### Changes that alter existing figures

Check these before re-running a figure made with a pre-release build.

- Dotplot and matrixplot draw the first category at the top, as in Scanpy.
  Use `swap_axes=True` rather than `coord_flip()`.
- Defaults follow Scanpy: rank-genes `n_genes=10`, rank matrixplot
  `standard_scale=None`, `plot_variance_ratio(n_pcs=30, log=False)`,
  `plot_highest_expr_genes` 30 genes.
- `plot_stacked_violin` fills each violin by group median expression and uses
  Scanpy's orientation.
- Without stored `adata.uns["<key>_colors"]`, discrete scales use Scanpy's
  default palette instead of plotnine hues.
- Inside `style_context`, `plot_features` and `plot_density` default to the
  style's sequential map (viridis) instead of magma.
- `plot_ma` renamed its column argument to `lfc_col`; `lfc` is now a
  fold-change threshold, as in `plot_volcano`.
- `plot_upset` sorts intersections largest first (`ascending=True` restores
  the old order).
- Box, violin and sina points are larger and more opaque by default.

### Scanpy spelling

- Every helper with `group_by=` also accepts Scanpy's `groupby=`; giving both
  with different values raises `TypeError`. The rank-genes helpers no longer
  warn for `groupby=`.
- Package warnings point at the caller's line, at any call depth.
- `plot_qc_violin` hides x tick labels with `group_by` on purpose; its
  docstring shows how to restore them.

### Differential expression, QC, pseudobulk and clustermap

- **Defaults now match scanpy.** `plot_rank_genes_dotplot` and
  `plot_rank_genes_matrixplot` use `n_genes=10`. `plot_rank_genes_matrixplot`
  uses `standard_scale=None` (raw group means). `plot_variance_ratio` uses
  `n_pcs=30, log=False`. `plot_highest_expr_genes` shows 30 genes by default.
- `plot_highest_expr_genes` takes `n_top` as an alias of `n`; giving both with
  different values raises. New `gene_symbols=` (label column from `var`) and
  `log=` (draws `log1p(percent)`, axis labelled so; zeros stay on the plot).
  The scaled-matrix warning now says "adata.raw is available: pass
  use_raw=True" when raw exists.
- `plot_rank_genes_dotplot` / `plot_rank_genes_matrixplot` gain `groups=`,
  `min_logfoldchange=` and `values_to_plot=` (`scores`, `logfoldchanges`,
  `pvals`, `pvals_adj`, `log10_pvals`, `log10_pvals_adj`). Selection happens at
  the DE-result boundary. `values_to_plot` draws the stored DE statistic, never
  an expression mean, and needs a result stored with `n_genes=adata.n_vars`.
  It cannot be combined with `standard_scale`.
- A gene that is a top marker of several groups now gets one column per marker
  group, labelled `GENE (group)`, with a `marker_group` column in `plot.data`.
  Values come from one computation per gene. Plotnine backend only.
- `n_genes <= 0`, non-integers and `None` raise `ValueError` in the rank-genes
  plots. Scanpy's bottom-ranked selection (negative `n_genes`) is unsupported.
- `plot_volcano` requires `group` to be one `str` (a list or `None` used to pool
  all groups). `color=`/`fill=` keywords are rejected because colour encodes
  direction. Default points are a little stronger (`size=1.5, alpha=0.8`).
- `plot_ma` signature changed: the column argument is now `lfc_col`
  (default `"log2FoldChange"`); `lfc` is a threshold (default `1.0`). A gene is
  significant when `padj < padj` and `|LFC| >= lfc`; points are coloured
  up/down/not significant like the volcano, with dashed lines at `+-lfc`.
  `plot.data` keeps a boolean `significant` column. Pass `lfc=0` for the old
  p-value-only rule. Rows with missing or non-positive `mean` are dropped with a
  warning giving the count; an all-dropped table raises. `data` must be a
  DataFrame; `label_top` must be a non-negative integer; `label` must exist.
- `plot_variance_ratio` validates `n_pcs` (positive integer).
- `plot_qc_violin` gains `log=` (`log1p(value)`, axis labelled) and rejects
  non-numeric metrics. `plot_qc_scatter` gains `alpha=` and `log=` (log10 axes;
  refuses non-positive values rather than dropping cells).
- `pseudobulk(use_raw=True)` now works when `adata.raw` has genes that
  `adata.var` lacks: counts and genes come from `adata.raw`, obs metadata is
  kept. A missing `adata.raw`, `sample_col` or `group_by` raises a clear error.
  Dropping profiles below `min_cells` now warns with the count.
- `plot_clustermap`: float32/int32 numeric annotations work. A clustering axis
  with fewer than two entries, or an empty matrix, raises `ValueError` before
  PyComplexHeatmap runs. Progress text is off (`verbose=0` unless passed). The
  fastcluster "Clustering large matrix" warning and the `set_bad` pending
  deprecation are filtered by message. The docstring says how to draw
  dendrograms (`row_dendrogram=True`, `col_dendrogram=True`).
- `plot_dendrogram` draws with `geom_path` in both orientations (the left
  orientation was a zig-zag). A single observed category, a non-categorical
  column and a missing column raise clear errors; unused categories are ignored.
- `plot_correlation` gains `use_rep=` / `n_pcs=` for scanpy's PCA-space
  statistic. The all-genes warning now appears only above 2000 genes.
- `plot_upset`: intersections are sorted largest first by default
  (`ascending=True` for the old order); all-empty sets and unknown
  `sort_subsets` raise. No set-size axis title is added (marsilea draws counts).
- Examples compute a fresh DE result into `rank_wilcoxon` and pass `key=`; the
  bundled `rank_genes_groups` is no longer overwritten.

### Distributions, markers and density

- `plot_violin`, `plot_box`, `plot_sina`, `plot_stacked_violin`, `plot_ridge`, `plot_expression_bar` accept numeric `adata.obs` columns in `genes` (bare name, `ggann.obs("x")` or `"obs:x"`). A bare name that is also a gene stays a gene.
- New `log=True` (violin, box, sina, stacked): plots `log1p(expression)`, zeros kept, axis labelled `log1p(expression)`. Values <= -1 raise.
- New `order=` (subset and order of groups) and `rotation=` on violin, box, sina, stacked violin. `plot_violin` gains `jitter_size`/`jitter_alpha`; `plot_stacked_violin` gains `stripplot`, `yticklabels`, `cmap`.
- `plot_stacked_violin` now fills each violin by group median expression (shared colour bar, `cmap="Blues"`), as scanpy does. Rows have no y tick labels by default. The group palette fill is gone.
- Grouped marker helpers: a gene under several labels is drawn once per label (no more ValueError); values come from one computation per gene. `plot_tracksplot` accepts a `{label: [genes]}` dict (tracks titled `gene [label]`). Non-dict `gene_groups` raises TypeError.
- `plot_tracksplot` draws `geom_rect` (about 12x to 23x faster, identical rectangles) and labels groups on the x axis.
- Seeded jitter: `random_state` (0 included) now seeds box/violin point jitter and sina spread; same positions on every draw/export. `None` is nondeterministic and no longer breaks a second draw.
- Default point size/alpha raised: box jitter 0.6/0.5, violin points 0.8/0.5, sina 0.8/0.6; sina legend keys enlarged.
- `stats=True` (violin, box) places the label 10% above the data and always adds headroom.
- `plot_sina`: `bins` must be a positive int; no fill-scale warning with `violin=False`.
- `plot_ridge`: validates `n_grid` (int >= 2), `scale` (finite > 0); new `bw_adjust`; raises above `MAX_RIDGE_GROUPS` (100) levels; warns for constant genes.
- `plot_expression_line(bins=N)` bins a numeric x. `plot_proportions` rejects `split_by == group_by`.
- `plot_density`: `downsample=` and `random_state=`; `joint=True` with one feature warns.

### Embeddings, matrix plots and the grammar

- `plot_embedding_density`: a single group no longer raises
  `Expected a 1D array`. Density is written by position, so interleaved groups and
  duplicate observation names stay aligned. Values are unchanged for every case
  that worked before. Points are drawn low to high density (scanpy order), so
  `plot.data` is sorted by density.
- `plot_embedding_density`: raises `ValueError` before any KDE when `group_by` has
  more than 100 levels (a numeric column is one level per value; 2,638 levels
  took 90 s).
- `plot_embedding_density`: new `groups=`, `background=` (grey context cells in
  every panel; layer data, not rows of `plot.data`, no part of any density),
  `components=`, `vmin=`, `vmax=`.
- `plot_embedding`: a missing category in `color=` keeps its cells, is drawn in
  `na_color`, and reads `NA` in the legend (was `nan`). Missing cells are drawn
  first. `plot.data` still holds the missing values.
- `plot_embedding`: new `color=[...]` (one panel per colour, `ncol=`),
  `components=`, `palette=`, `color_map=`, `vmin/vmax/vcenter` (numbers or
  `"p99"`; values clip), `groups=`, `sort_order=`, `na_color=`, `title=`,
  `legend_loc=` (`"right margin"`, `"right"`, `"on data"`, `"none"`).
  `backend="matplotlib"` rejects these with a `ValueError`.
- `plot_embedding(label=True, split_by=<same column>)` no longer raises
  `cannot insert ... already exists`, and labels land in their own facet.
- `plot_embedding`: new `method="binned"` for the no-colour point density
  (histogram2d + Gaussian filter, sampled at each cell). Opt-in; the legend reads
  `density (binned)`. The exact KDE stays the default and warns above 20,000
  cells, suggesting `method="binned"`, `pointdensity=False` or `downsample=`.
  The statsmodels `FutureWarning` about entropy-initialized generators is gone.
- `plot_features`: with an obs column among the features (`color_scale="auto"`),
  each feature gets its own `plot_embedding` panel and colour bar, composed into
  a grid. Genes only keep the single shared scale (`color_scale="shared"`
  forces it, `"free"` forces panels). Values are never rescaled. New
  `components=`, `vmin/vmax/vcenter`, `sort_order=`, `na_color=`. The return
  value is a composition in the free case.
- `plot_dotplot` / `plot_matrixplot`: the first category is drawn at the top, as in
  scanpy, in both backends. Data and category order are unchanged. Use
  `swap_axes=True` rather than `coord_flip()` (a flip reverses the displayed
  order).
- `plot_dotplot` / `plot_matrixplot` / `plot_heatmap`: new `swap_axes`,
  `vmin`/`vmax`, and dict `genes` (`{marker group: [genes]}`: one panel per
  group; a repeated gene appears once per occurrence, labelled by group, with
  values from one computation; plotnine backend only). `plot_dotplot` also gets
  `dot_min`/`dot_max`. `cmap` now defaults to `None` and resolves to the old
  default; with `standard_scale="zscore"` it resolves to `RdBu_r` with
  zero-centred limits (was a sequential map).
- `plot_dotplot` / `plot_matrixplot`: a `UserWarning` when `group_by` has more than
  100 levels (a numeric column is one level per value).
- `plot_heatmap` / `plot_tracksplot`: cells keep their input order within a group
  (stable sort in the shared `_cell_rank`).
- `gganndata`: fields named later with `+ aes(...)`, `+ geom_*(aes(...))` or
  `+ facet_*(...)` are resolved from the AnnData within the original
  `max_matrix_values` budget. Expressions such as `np.log1p(n_counts)`,
  `factor(louvain)` and `Q("HLA-DRA")` project the fields they use. Any column of
  any `obsm` array resolves (`PC_5`, DataFrame-backed `obsm` by column position).
  `add_theme=False` keeps plotnine's default `theme_gray`; the return type is the
  same (`GGAnnPlot`).

### Themes, palettes and publication export

- `set_theme`, `theme_ggann` and `sizes.update` reject a `base_size` that is not a finite positive number
  (`ValueError`) before any global state changes.
- `reset_theme()` now also resets `ggann.sizes` to a base size of 11.
- A theme or style font family that is not installed (for example `set_theme(family="NoSuchFont")` or
  Arial in a fresh Linux install) now emits one `UserWarning` per family and process, then falls back as before.
- `scale_color_expression` and `scale_fill_expression` take `cmap=None` by default. `None` means Reds
  (colour) or viridis (fill) outside a style context and the style's `sequential_cmap` inside one. An explicit
  name always wins, including `"Reds"` inside `style_context` (it used to be replaced by viridis).
- `plot_features`, `plot_dotplot`, `plot_matrixplot`, `plot_embedding_density`, `plot_heatmap` and `plot_density`
  take `cmap=None` with the same rule (family defaults stay magma, Reds, viridis, viridis, viridis, magma).
  Inside `style_context`, `plot_features` and `plot_density` now default to the style's sequential map
  (viridis) instead of magma. An explicit `cmap` with `standard_scale="zscore"` in a style context now
  uses that map with zero-centred limits instead of the diverging style palette.
  Precedence: explicit `cmap` > style `sequential_cmap` > family default; `plot_embedding`'s `low`/`high` only
  drive the legacy two-colour gradient and are ignored in a style context.
- `scale_color_obs` / `scale_fill_obs`: numeric observation columns get a continuous colour-map scale
  (viridis, or the style's sequential map); boolean and categorical columns stay discrete; a missing column raises
  `KeyError` at construction. `obs_colors` still returns a dict or `None`.
- With no stored `adata.uns["<col>_colors"]` and no publication style, discrete observation scales and the
  Matplotlib backend now use scanpy's default palette (20, 28 or 102 colours) instead of plotnine hues, so
  colours match `sc.pl.umap`. More than 102 categories keep hues. `adata.uns` is never modified.
- `plot_clustermap` annotation bars reuse stored (or publication, or scanpy default) category colours when
  every annotation column is a discrete `obs` column.
- Publication qualitative palettes beyond the style vocabulary are extended from a fixed colour sequence
  (the old resampled `turbo` palette is gone), so 20, 21 and 22 categories never recolour earlier ones.
- One `UserWarning` per palette when stored or generated publication colours fall below 5 CIEDE2000 under a
  colour-vision simulation (only the first eight colours are a measured guarantee), and when stored colours are
  duplicated. Invalid stored colours now warn and fall back to the default palette in the non-publication path too.
- `save_publication`:
  - accepts any positive integer `dpi` (150 works); bool, zero, negative, non-finite and non-integer values raise.
    `PublicationStyle(dpi=...)` still accepts only 300 or 600. Raster size is `floor(size_in * dpi)` pixels.
  - warns when text, a legend or a colourbar extends beyond the final canvas or is clipped by its axes. The canvas
    is never cropped or enlarged.
  - treats only svg/pdf/png/tiff/tif as a suffix. `fig1.v2` writes `fig1.v2.svg/.pdf/.png`. A suffix that
    is missing from an explicit `formats` raises (`x.png` with `formats=("svg",)` used to write `x.svg`).
    `.jpg`, `.eps` and similar still raise. `height` stays required.
- `compose(...).save(..., dpi=150)` works (any positive integer, default 300).
- `compose` and `tag_panels` raise `TypeError` for non-plot panels.
- The Matplotlib dotplot/matrixplot left margin is measured from the rendered y labels (and re-measured at export
  size) instead of a character-count estimate.
