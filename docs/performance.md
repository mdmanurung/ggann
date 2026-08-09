# Performance and Scanpy comparison

```{warning}
The frozen release candidate fails the measured Scanpy speed and memory gates.
The results below do **not** support a "blazing fast" claim. Publication remains
blocked until a repeated full matrix passes every gate.
```

## Current evidence status

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

### Frozen-candidate regression check

The exact pre-optimization source (`b8c9eb…`) was rerun with the same fixture,
environment, seed, and seven repetitions. The candidate source is `7e059e8…`.

| Workload | Preparation change | Construction change | Render/save change | End-to-end change |
|---|---:|---:|---:|---:|
| Embedding | -19.9% | -36.7% | -29.2% | -18.1% |
| Dotplot | -2.5% | -13.2% | -1.5% | -1.0% |
| Matrixplot | -12.7% | -20.2% | +3.0% | -1.5% |

All primary end-to-end and preparation timings pass the 5% regression limit.
A single fresh-child matrixplot end-to-end RSS probe increased from 16.19 to
17.16 MiB (+6.0%), while profiler allocations were essentially unchanged. This
memory observation is reported as an unresolved baseline-regression blocker,
not dismissed as a pass.

### Mandatory release gates

| Primary workload | Preparation gate (≥2×) | E2E gate (≥0.91×) | Peak RSS gate (≤1.10×) | Gate |
|---|---:|---:|---:|---|
| Embedding | 0.33× — fail | 0.55× — fail | 5.19× — fail | **FAIL** |
| Dotplot | 0.35× — fail | 0.50× — fail | 1.40× — fail | **FAIL** |
| Matrixplot | 0.23× — fail | 0.53× — fail | 2.22× — fail | **FAIL** |
| Primary-suite geometric mean | — | 0.526× — fail | — | **FAIL** |

The performance verdict remains `NOT READY`; these results do not support a
"blazing fast" claim.

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

The smoke artifacts are disposable CI checks; they are not release evidence.

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

The full dense/CSR/CSC × `.X`/layer/raw matrix remains required before a release
can pass, but it cannot overturn failures already present in the primary CSR
suite.

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
  --label candidate \
  --output benchmarks/results/candidate.json
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
  --label candidate-large \
  --output benchmarks/results/candidate-large.json
```

Compare the same environment and fixture before and after optimization:

```bash
python benchmarks/compare_results.py \
  benchmarks/results/baseline.json \
  benchmarks/results/candidate.json \
  --fail-regression-pct 5 \
  --fail-on-output-change \
  --output benchmarks/results/baseline-v-candidate.md
```

These commands measure ggann regressions. A release also requires the matched
Scanpy comparison artifact described above; a ggann-only result cannot satisfy
the Scanpy gates.

## Correctness before speed

Prepared summaries must match within documented numeric tolerances, including
group order, gene order, means, fractions expressing, and the represented cell
set. After optimization, compare representative rendered artifacts and prepared
tables against the baseline:

```bash
PYTHONPATH=/path/to/baseline python benchmarks/check_invariance.py snapshot \
  --output /tmp/ggann-baseline
PYTHONPATH=src python benchmarks/check_invariance.py snapshot \
  --output /tmp/ggann-candidate
python benchmarks/check_invariance.py compare \
  /tmp/ggann-baseline /tmp/ggann-candidate
```

The `pbmc68k_reduced` artifact check may use a pre-populated Scanpy cache for
offline execution. Deterministic synthetic fixtures remain the network-free
fallback when that dataset is unavailable.

## Release gates

The release passes only when all of the following are true:

- primary large sparse embedding, dotplot, and matrixplot preparation is at
  least 2× faster than Scanpy by median time;
- no primary end-to-end workload is more than 10% slower than Scanpy;
- the end-to-end geometric mean is faster across the primary suite;
- large sparse peak memory is no more than 10% worse than Scanpy;
- no ggann workload regresses by more than 5% from the frozen ggann baseline
  without a documented correctness justification.

Store the raw JSON, comparison report, invariance report, exact commands, and
environment provenance together. The final values and artifact paths belong in
`RELEASE_READINESS.md` as well as the status table at the top of this page.
