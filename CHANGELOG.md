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
