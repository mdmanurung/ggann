# ggann 0.1.0 verification record

Assessment date: 2026-08-09

This record captures the local evidence for ggann's first public release. The
verified scope is the AnnData grammar, all 30 plotting helpers, the optional
publication design system, exact multi-format export, executable documentation,
and source/distribution reproducibility.

## Outcome

| Area | Result | Evidence |
|---|:---:|---|
| Python API and helpers | **PASS** | 386 tests on each of Python 3.12, 3.13, and 3.14; all 28 plotnine helpers and both grid-backend exceptions exercised |
| Static quality | **PASS** | Ruff check and format verification clean; Pyright reports zero errors, warnings, or information messages for `src/ggann` |
| Documentation | **PASS** | Offline Sphinx build succeeds with warnings as errors and executes all six workflows |
| Scientific invariance | **PASS** | Frozen `30e3359` prepared tables match at `rtol=1e-6`, `atol=1e-7`; mappings, labels, scales, populations, category order, and AnnData fingerprints match exactly |
| Exploratory rendering | **PASS** | Representative legacy embedding, dotplot, matrixplot, and heatmap PNGs are pixel-identical; default SVG output remains vector |
| Publication rendering | **PASS** | Exact canvases, final-size text/overlap checks, vector/raster layering, fonts, palette accessibility, object cleanup, and real PBMC output contracts pass |
| Performance regression | **PASS** | 18 timing/RSS checks against `30e3359` pass the 5% limit with seven repeated samples and seven isolated memory children |
| Distribution | **PASS** | Offline wheel and sdist build; strict Twine checks; wheel import outside the checkout; `pip check`; source/wheel output identity |

## Runtime matrix

All runs used the official annplyr 0.3.0 and plotnine-extra 0.3.1 source
snapshots, plotnine 0.15.7, Scanpy 1.12.3, AnnData 0.13.2, NumPy 2.4.6,
pandas 3.0.5, SciPy 1.18.0, and Matplotlib 3.11.1.

| Python | Command | Result |
|---|---|---:|
| 3.12.13 | `PYTHONPATH=src python -m pytest -q` | 386 passed |
| 3.13.11 | `PYTHONPATH=src python -m pytest -q` | 386 passed |
| 3.14.5 | `PYTHONPATH=src python -m pytest -q` | 386 passed |

Expected third-party warnings cover PyComplexHeatmap's pending colormap
deprecation, deliberately duplicated observation names, tiny density groups,
and fixture-specific scaled-expression guidance. No test is skipped or marked
expected-failure.

## Scientific and visual invariance

The representative real-data audit uses
`scanpy.datasets.pbmc68k_reduced`. It compares frozen and candidate:

- exact AnnData SHA-256 fingerprints before and after construction/rendering;
- resolved and aggregated tables, including structure and categorical order;
- plot-ready embedding, dotplot, matrixplot, and heatmap tables;
- mappings, labels, scales, visible text, represented point/tile counts, and
  rasterization state;
- PNG channel differences and SVG raster-image presence.

Every prepared result passes the numeric tolerances. All four PNG pairs have
maximum and mean channel differences of zero and infinite PSNR. The complete
report is [publication-legacy-invariance.md](benchmarks/results/publication-legacy-invariance.md).

Publication tests additionally cover 89 × 70 mm and 183 × 120 mm canvases,
one-pixel raster dimension tolerance, 0.5-point canvas containment, less than
5% unintended text overlap, a 5-point minimum visible font, exact represented
populations, selective layer rasterization, editable SVG text, embedded PDF
TrueType/CID Type 2 fonts, 20 repeated exports with less than 1 MiB retained
Python allocation growth, and exception-path cleanup.

The core eight-colour vocabulary passes minimum CIEDE2000 separation of 10 in
normal vision and 5 under protanopia, deuteranopia, and tritanopia simulations.
Sequential luminance and both diverging arms are monotonic.

## Real PBMC evidence

The claim-driven workflow represents all 700 cells in the bundled PBMC68k
subset and combines lineage UMAP, five-gene dotplot, CD3D distribution, and
lineage-by-phase composition panels. It writes:

- exploratory and publication PNGs;
- editable publication SVG and PDF;
- exact-size 2161 × 1417 PNG and TIFF for 183 × 120 mm at 300 DPI;
- an accessibility report;
- a manifest containing the claim, panel map, `n`, summary definitions, style,
  palette, software versions, input fingerprint, canvas, DPI, and hashes.

The retained local evidence is under `examples/_output/pbmc-publication/`; the
two review images are versioned under `docs/images/`. The workflow itself is
`examples/vignettes/05_publication_panels.py` and can recreate every artifact
without network access.

## Performance regression

Frozen commit `30e3359` and the candidate were run sequentially on the same
extended 20,000 × 10,000 CSR fixture with one BLAS/OpenMP thread, identical
seeds/caches, seven repeated timings, and seven fresh-process RSS measurements
per library and stage.

| Workload | Preparation | Construction | Render | End to end | Peak RSS | Retained RSS |
|---|---:|---:|---:|---:|---:|---:|
| Embedding | 0.963× | 0.954× | 0.998× | 0.993× | 1.002× | 0.960× |
| Dotplot | 1.003× | 1.002× | 1.000× | 0.996× | 0.932× | 0.932× |
| Matrixplot | 0.993× | 1.009× | 0.993× | 1.000× | 0.934× | 0.934× |

Ratios are candidate divided by frozen; values at or below 1.05 pass. Raw
samples and the machine-readable comparison are versioned as
`benchmarks/results/publication-*.json`, with the concise table in
[publication-regression.md](benchmarks/results/publication-regression.md).

The separate matched Scanpy benchmark favours Scanpy for speed and memory on
the measured large sparse workloads. The documentation therefore makes no
blanket speed claim and positions ggann around explicit source selection,
plotnine composability, and exact publication export.

## Documentation workflows

The warning-as-error offline build executes six scenario-led examples:

1. a real PBMC marker-review figure from a familiar Scanpy analysis;
2. a condition-associated MKI67 question expressed through the grammar;
3. sparse and read-only backed marker review;
4. a condition-aware custom annplyr summary beside a standard helper;
5. a real PBMC claim-to-publication workflow;
6. an evidence-linked Scanpy trade-off comparison.

The first and fifth use bundled real PBMC data. The focused workflows use a
deterministic PBMC-like control/stimulation fixture so documentation remains
fast and network-free.

## Distribution verification

`uv build --offline` produced `ggann-0.1.0-py3-none-any.whl` and
`ggann-0.1.0.tar.gz`.

`twine check --strict` passes both. The wheel was installed over the package in
a dependency-complete environment copied to `/tmp/ggann-wheel-verify`, then
imported from that environment's `site-packages` while the working directory
was `/tmp`. `pip check` reports no broken requirements.

Wheel-installed legacy PNGs are pixel-identical to source, with identical
prepared tables, mappings, labels, scales, artists, and SVG vector contracts.
The wheel-installed real PBMC exploratory PNG, publication PNG, and TIFF are
also pixel-identical to source. Publication manifests match, SVG text remains
editable with no image layer, and both PDFs reject Type 3 fonts. See
[publication-wheel-invariance.md](benchmarks/results/publication-wheel-invariance.md).

The build used cached declared build requirements because the verification
environment had no network access. Public-index dependency resolution is an
operational upload-time check and was not inferred from the offline build.

## Reproduction commands

```bash
python scripts/run_pyright.py
ruff check src tests benchmarks docs/extensions examples scripts
ruff format --check src tests benchmarks docs/extensions examples scripts
PYTHONPATH=src MPLBACKEND=Agg python -m pytest -q
GGANN_DOCS_OFFLINE=1 MPLBACKEND=Agg \
  sphinx-build -W --keep-going -b html docs docs/_build/html
uv build --offline --out-dir dist
twine check --strict dist/*
```

The exact benchmark and invariance commands are documented in
`docs/performance.md` and retained in the JSON metadata.
