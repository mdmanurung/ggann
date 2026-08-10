# Performance and Scanpy comparison

```{warning}
The matched benchmark below favours Scanpy for speed and memory on the measured
large sparse workloads. Choose ggann for its tidy AnnData-to-plotnine workflow,
composable grammar, and exact publication export—not for a blanket speed claim.
```

## Evidence and interpretation

### What the refined profile found

`benchmarks/profile_hotspots.py` records cProfile calls, tracemalloc lines,
sampled peak/retained RSS, DataFrame sizes, sparse-column counts, and every
public annplyr `to_df` request. The same 20,000 × 10,000 CSR fixture was profiled
before and after the changes; these diagnostic medians use three uninstrumented
warm calls, while cProfile and tracemalloc run separately.

| Workload and stage | Before | After | Change | annplyr calls |
|---|---:|---:|---:|---:|
| Embedding preparation | 20.3 ms | 14.5 ms | -28.5% | 2 → 1 |
| Embedding public construction | 27.9 ms | 15.6 ms | -43.9% | 2 → 1 |
| Embedding `ggplot.draw()` | 365.1 ms | 346.0 ms | -5.2% | 0 → 0 |
| Embedding PNG save | 589.4 ms | 399.9 ms | -32.1% | 0 → 0 |
| Dotplot preparation | 56.2 ms | 52.3 ms | -6.8% | 2 → 1 |
| Dotplot end to end | 758.1 ms | 734.3 ms | -3.1% | 2 → 1 |
| Matrixplot preparation | 65.4 ms | 53.5 ms | -18.1% | 2 → 1 |
| Matrixplot end to end | 685.1 ms | 636.6 ms | -7.1% | 2 → 1 |

The concrete profile attributes the remaining preparation cost to annplyr's
public projected extraction: sparse source adaptation, variable-name position
resolution over 10,000 names, projected sparse-frame construction, and selector
evaluation. ggann's sparse-native group reduction is no longer the dominant
path. Matrixplot's previous `DataFrame.melt` cost was removed by constructing
the compact group-by-gene long frame directly.

Rendering remains the larger limit. Embedding spends most draw time in
`geom_point.draw_unit`, Matplotlib `Axes.scatter`, colour conversion, and
collection construction. Dotplot and matrixplot are dominated by plotnine
scale/guide/layout work, theme application, tick creation, font lookup, and
Matplotlib text/artist drawing. These are costs of retaining an ordinary,
composable plotnine object.

Implemented optimizations are deliberately narrow:

- compatible expression/`obs` and `obsm`/colour/facet projections cross the
  annplyr boundary once;
- cumulative budgets are still checked before the first read, typed errors
  propagate, and output is split/assembled by position;
- group-by-gene data stays wide until the one required compact reshape;
- plot components are composed in one plotnine addition, avoiding repeated
  deep copies of a fresh plot and theme;
- embedding points use a zero-width redundant same-colour outline with a size
  offset that preserves their rendered diameter.

`geom_raster` produced pixel-identical matrixplot PNGs and reduced draw time in
an experiment, but it was not adopted because it would change vector-export
semantics and did not improve PNG save time. No cells, genes, legends, facets,
statistics, resolution, or output format were removed or downsampled.

### Repeated large-CSR primary result

The versioned JSON retains all seven timing samples, cold calls, fresh-child
memory probes, complete provenance, payload tolerances, and input fingerprints:
`benchmarks/results/scanpy-extended-csr.json`.

```{include} _includes/scanpy-extended-csr.md
```

Fresh-process imports took 3.63 seconds for ggann and 2.97 seconds for Scanpy.
The run used Python 3.12.13 on a 24-logical-core Intel Xeon E5-2690 v3 with one
BLAS/OpenMP thread, Scanpy 1.12.3, annplyr 0.3.0, AnnData 0.13.2, NumPy 2.4.6,
pandas 3.0.5, SciPy 1.18.0, plotnine 0.15.7, plotnine-extra 0.3.1, and Matplotlib
3.11.1.

### Frozen first-release self-regression check

The publication implementation was compared with frozen commit `30e3359` using
the same extended CSR fixture, environment, seed, and seven repeated samples.
Peak and retained RSS are medians from seven independent fresh children per
library and stage.

| Workload | Preparation | Construction | Render | End to end | Peak RSS | Retained RSS |
|---|---:|---:|---:|---:|---:|---:|
| Embedding | -3.7% | -4.6% | -0.2% | -0.7% | +0.2% | -4.0% |
| Dotplot | +0.3% | +0.2% | +0.0% | -0.4% | -6.8% | -6.8% |
| Matrixplot | -0.7% | +0.9% | -0.7% | +0.0% | -6.6% | -6.6% |

All 18 checks pass the 5% regression limit. The largest timing ratio is 1.009×
for matrixplot construction; the largest memory ratio is 1.002× for embedding
peak RSS. The versioned `publication-baseline-30e3359.json`,
`publication-candidate.json`, and `publication-regression.json` files under
`benchmarks/results/` retain every sample and source-tree digest.

### Conservative comparison thresholds

| Primary workload | Preparation gate (≥2×) | E2E gate (≥0.91×) | Peak RSS gate (≤1.10×) | Gate |
|---|---:|---:|---:|---|
| Embedding | 0.33× — fail | 0.55× — fail | 5.19× — fail | **FAIL** |
| Dotplot | 0.35× — fail | 0.50× — fail | 1.40× — fail | **FAIL** |
| Matrixplot | 0.23× — fail | 0.53× — fail | 2.22× — fail | **FAIL** |
| Primary-suite geometric mean | — | 0.526× — fail | — | **FAIL** |

These thresholds are deliberately demanding and are not met by this run. The
result is a limitation on performance claims, not on the correctness or
ergonomics of the plotting API.

## Matched plot pairs

The comparison suite must run the nearest semantic equivalents on the same
immutable `AnnData`:

| ggann | Scanpy |
|---|---|
| `plot_embedding` | `scanpy.pl.embedding` |
| `plot_dotplot` | `scanpy.pl.dotplot` |
| `plot_matrixplot` | `scanpy.pl.matrixplot` |
| `plot_violin` | `scanpy.pl.violin` |
| `plot_stacked_violin` | `scanpy.pl.stacked_violin` |
| `plot_tracksplot` | `scanpy.pl.tracksplot` |
| `plot_highest_expr_genes` | `scanpy.pl.highest_expr_genes` |
| `plot_rank_genes_dotplot` | `scanpy.pl.rank_genes_groups_dotplot` |
| `plot_rank_genes_matrixplot` | `scanpy.pl.rank_genes_groups_matrixplot` |

A pair is excluded rather than forced when it cannot represent the same cells,
variables, grouping, statistic, and output quality.

## Measure stages separately

Every accepted pair records five non-overlapping views of cost:

1. **Preparation**: source resolution, projection, joins, aggregation, and the
   final plot-ready table.
2. **Construction**: creation of the plotting object from already prepared
   data.
3. **Rendering**: drawing an already constructed object with a warmed renderer.
4. **End to end**: preparation, construction, and an equivalent headless save.
5. **Memory**: sampled peak RSS and retained RSS after the output is released
   and garbage collected.

Preparation speedups and renderer speedups must be reported separately.
Object construction alone is not an end-to-end result.

## Fixtures and controls

Required cases cover dense, CSR, and CSC data at small, medium, and realistic
large scales. The large sparse primary cases use at least 20,000 observations,
10,000 variables, 32 requested genes, and 32 groups. Cases include categorical
and continuous embedding colour, several genes and groups, `.X`, `.raw`, and a
named layer.

For every pair:

- reuse the exact same immutable fixture and verify its fingerprint;
- pin the random seed and BLAS/OpenMP thread counts;
- use Matplotlib's headless backend and the same width, height, DPI, and output
  format;
- warm font and renderer caches before repeated samples;
- measure cold calls and peak memory in fresh child processes;
- retain every timing sample, not only its summary;
- report median, range or confidence interval, speedup ratio, peak RSS, and
  retained RSS;
- record OS, CPU, Python, ggann, Scanpy, annplyr, AnnData, NumPy, pandas, SciPy,
  plotnine, and Matplotlib versions.

The result is invalid if either library plots fewer cells, silently downsamples,
changes the requested statistic, densifies unrelated sparse columns, or writes a
lower-quality figure.

## Reproduce the benchmarks

### Matched Scanpy smoke

First reproduce the matched smoke comparison:

```bash
NUMBA_CACHE_DIR=/tmp/ggann-numba-cache \
MPLCONFIGDIR=/tmp/ggann-mpl-cache \
PYTHONPATH=src python benchmarks/compare_scanpy.py \
  --preset smoke \
  --formats csr \
  --workloads primary \
  --sources x \
  --repeats 2 \
  --seed 20260809 \
  --output /tmp/ggann-vs-scanpy-smoke.json \
  --report /tmp/ggann-vs-scanpy-smoke.md
```

The smoke artifacts are disposable CI checks; they are not evidence for a
general performance claim.

Reproduce the committed primary result on a quiet dedicated machine:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
NUMBA_CACHE_DIR=/tmp/ggann-benchmark-numba \
MPLCONFIGDIR=/tmp/ggann-benchmark-mpl MPLBACKEND=Agg PYTHONHASHSEED=0 \
PYTHONPATH=src python benchmarks/compare_scanpy.py \
  --preset extended \
  --formats csr \
  --workloads primary \
  --sources x \
  --repeats 7 \
  --seed 20260809 \
  --include-cold-start \
  --isolated-memory-stages preparation,end_to_end \
  --output benchmarks/results/scanpy-extended-csr.json \
  --report benchmarks/results/scanpy-extended-csr.md

python benchmarks/render_scanpy_vignette.py \
  benchmarks/results/scanpy-extended-csr.json \
  docs/_includes/scanpy-extended-csr.md
```

Run the concrete hotspot profiler separately; its instrumented timings are
diagnostic and never substituted for the matched benchmark:

```bash
PYTHONPATH=src python benchmarks/profile_hotspots.py \
  --preset extended --format csr --repeats 3 \
  --git-revision "$(git rev-parse HEAD)" \
  --output benchmarks/results/profile-hotspots-after.json \
  --profile-dir /tmp/ggann-profile-pstats
```

The full dense/CSR/CSC × `.X`/layer/raw matrix is required before making a broad
cross-library performance claim, but it cannot overturn failures already
present in the primary CSR suite.

### ggann baseline and regression harness

The repository's ggann-only harness is used separately to detect regressions.

It runs each case in a fresh child process and keeps raw samples, output
fingerprints, stage sizes, dependency provenance, and RSS measurements:

```bash
PYTHONPATH=src python benchmarks/run_benchmarks.py \
  --preset standard \
  --formats dense,csr,csc \
  --workloads all \
  --include-rendering \
  --repeats 5 \
  --seed 20260808 \
  --label current \
  --output benchmarks/results/current.json
```

Run the realistic sparse suite separately:

```bash
PYTHONPATH=src python benchmarks/run_benchmarks.py \
  --preset extended \
  --formats csr,csc \
  --workloads all \
  --include-rendering \
  --repeats 5 \
  --seed 20260808 \
  --label current-large \
  --output benchmarks/results/current-large.json
```

Compare the same environment and fixture before and after optimization:

```bash
python benchmarks/compare_results.py \
  benchmarks/results/baseline.json \
  benchmarks/results/current.json \
  --fail-regression-pct 5 \
  --fail-on-output-change \
  --output benchmarks/results/baseline-v-current.md
```

These commands measure ggann self-regressions. A cross-library performance
claim also requires the matched Scanpy artifact described above; a ggann-only
result cannot establish comparative speed.

## Correctness before speed

Prepared summaries must match within documented numeric tolerances, including
group order, gene order, means, fractions expressing, and the represented cell
set. After optimization, compare representative rendered artifacts and prepared
tables against the baseline:

```bash
PYTHONPATH=/path/to/baseline python benchmarks/check_invariance.py snapshot \
  --output /tmp/ggann-baseline
PYTHONPATH=src python benchmarks/check_invariance.py snapshot \
  --output /tmp/ggann-current
python benchmarks/check_invariance.py compare \
  /tmp/ggann-baseline /tmp/ggann-current
```

The benchmarks keep the `pbmc68k_reduced` artifact check so recorded timings
stay comparable against the published results, and because Scanpy's bundled copy
needs no network. The documentation, vignettes, and example scripts use
`pbmc3k_processed`, which `scripts/fetch_datasets.py` downloads once into
`data/`; with `GGANN_DOCS_OFFLINE=1` a missing cache raises rather than
downloading.

## Performance acceptance criteria

ggann uses the following criteria before making broad speed or memory claims:

- primary large sparse embedding, dotplot, and matrixplot preparation is at
  least 2× faster than Scanpy by median time;
- no primary end-to-end workload is more than 10% slower than Scanpy;
- the end-to-end geometric mean is faster across the primary suite;
- large sparse peak memory is no more than 10% worse than Scanpy;
- no ggann workload regresses by more than 5% from the frozen ggann baseline
  without a documented correctness justification.

Store the raw JSON, comparison report, invariance report, exact commands, and
environment provenance together so every published number remains auditable.
