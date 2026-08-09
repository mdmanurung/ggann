# Changelog

All notable changes to ggann are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and releases use
[semantic versioning](https://semver.org/).

## [Unreleased]

### Added

- A plotnine-native grammar for resolving AnnData observation fields, genes,
  layers, raw expression, and embeddings.
- High-level helpers covering embeddings, expression summaries,
  distributions, differential expression, composition, quality control, and
  figure composition.
- Reproducible performance and downstream-artifact invariance tooling.
- Executable, offline documentation vignettes and explicit public API
  contracts.
- A concrete hotspot profiler and an executable, evidence-linked comparison of
  ggann + annplyr with standard Scanpy plotting.
- Optional cumulative `max_matrix_values` materialization budgets on
  `gganndata()` and `embedding_coords()`.

### Changed

- Runtime dependencies now use release version ranges instead of Git commit
  references.
- AnnData extraction targets annplyr 0.3's public projected-read interface.
- Grouped summaries now use one projected annplyr read followed by a
  sparse-native vectorized reduction; missing grouping keys are omitted and
  expression `NaN` values follow pandas/Scanpy mean and detection semantics.
- Float64 group accumulation improves numerical stability while preserving
  downstream prepared data within the documented `rtol=1e-6`, `atol=1e-7`
  tolerance.
- Compatible coordinates, colour/facet metadata, expression, and grouping
  fields are fused into fewer public annplyr projections while retaining
  cumulative preflight budgets and positional reconstruction.
- Primary plotnine helpers avoid repeated plot/theme deep copies; embedding
  points avoid redundant same-colour outline rendering without changing their
  displayed diameter or return type.

## [0.1.0] - Unreleased

Initial public release candidate.

[Unreleased]: https://github.com/mdmanurung/ggann/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/mdmanurung/ggann/releases/tag/v0.1.0
