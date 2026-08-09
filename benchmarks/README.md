# Benchmarks

The benchmark runner measures AnnData extraction, aggregation, plot preparation,
and optional rendering. It generates seeded synthetic data and runs every case in
a fresh child process.

## Compare with Scanpy

`compare_scanpy.py` runs matched ggann and Scanpy plots on the same seeded,
fingerprinted AnnData fixture. It validates the prepared values before admitting
any speed ratio and checks that plotting did not change the input fingerprint.
The primary smoke comparison is:

```bash
NUMBA_CACHE_DIR=/tmp/ggann-numba-cache \
MPLCONFIGDIR=/tmp/ggann-mpl-cache \
PYTHONPATH=src python benchmarks/compare_scanpy.py \
  --preset smoke \
  --formats csr \
  --workloads primary \
  --sources x \
  --repeats 3 \
  --output benchmarks/results/scanpy-smoke.json \
  --report benchmarks/results/scanpy-smoke.md
```

`primary` covers categorical embedding, dotplot, and matrixplot. `all` adds
continuous and gene-coloured embeddings, violin, stacked violin, tracksplot,
highest-expressed genes, and the fair ranked-gene dotplot/matrixplot pairs. Use
`--list-workloads` to print the individual names. Expression sources can be
`x`, `layer`, `raw`, or `all`; named-layer cases use the deterministic `counts`
fixture layer. Scanpy's highest-expression helper has no `use_raw` argument, so
that unmatched raw case is explicitly skipped and recorded rather than adapted
to a different AnnData object.

Every case records raw cold-call and repeated samples for four stages:

- `preparation`: native extraction/statistical preparation, checked numerically;
- `construction`: the public plotting call with display disabled;
- `render`: matched-size/DPI PNG save from a materialized figure;
- `end_to_end`: public plotting call through PNG save.

Only `preparation` and `end_to_end` rows with `comparability=pass` may support a
cross-library speed claim. Construction is diagnostic because Scanpy eagerly
creates matplotlib artists for several functions while ggann returns a deferred
plotnine object. Rendering uses the same 6 by 4.5 inch, 80 DPI, Agg-backed PNG
target. Raw samples include sampled peak and retained RSS; short-lived
allocations below the sampling interval can still be missed.

The release-sized gate command is intentionally opt-in:

```bash
NUMBA_CACHE_DIR=/tmp/ggann-numba-cache \
MPLCONFIGDIR=/tmp/ggann-mpl-cache \
PYTHONPATH=src python benchmarks/compare_scanpy.py \
  --preset extended \
  --formats csr \
  --workloads primary \
  --sources x \
  --repeats 7 \
  --include-cold-start \
  --isolated-memory-stages preparation,end_to_end \
  --output benchmarks/results/scanpy-extended-csr.json \
  --report benchmarks/results/scanpy-extended-csr.md
```

The isolated-memory option creates one fresh child per library and selected
stage. Imports and fixture construction precede the RSS baseline, and the
measured call is the child's first preparation or plotting call. This avoids
allocator reuse from earlier stages obscuring peak RSS.

The JSON evaluates the cross-library preparation, end-to-end geometric-mean,
per-case slowdown, and peak-memory gates only when every extended sparse primary
case is present and comparable. The separate before/after ggann regression
runner below remains authoritative for the 5% self-regression gate; neither
runner alone is a complete performance-release verdict.

The repository retains the accepted raw primary samples in
`benchmarks/results/scanpy-extended-csr.json` and the exact frozen candidate in
`scanpy-extended-csr-before.json`. Regenerate the documentation table rather
than copying medians by hand:

```bash
python benchmarks/render_scanpy_vignette.py \
  benchmarks/results/scanpy-extended-csr.json \
  docs/_includes/scanpy-extended-csr.md
```

`cold` is the first workload call after imports and fixture construction. Repeated
timings reuse the same immutable AnnData object. They therefore measure repeated
package use without including data generation or import time.

## Profile concrete hotspots

`profile_hotspots.py` uses the same deterministic fixture but serves diagnosis,
not cross-library claims. It records unprofiled warm samples, concrete cProfile
functions, tracemalloc peaks and top lines, sampled RSS, DataFrame sizes, sparse
column counts, and every public annplyr extraction call. Raw pstats files are
written separately so they can be opened with `python -m pstats` or a profile
viewer:

```bash
PYTHONPATH=src python benchmarks/profile_hotspots.py \
  --preset extended --format csr --repeats 3 \
  --git-revision "$(git rev-parse HEAD)" \
  --output benchmarks/results/profile-hotspots-after.json \
  --profile-dir /tmp/ggann-profile-pstats
```

The versioned before/after JSON documents retain samples and top-function data.
cProfile and tracemalloc timings are instrumented separately and must not be
reported as ordinary benchmark timings.

## Run

Activate an environment with ggann and its development dependencies, then run:

```bash
PYTHONPATH=src python benchmarks/run_benchmarks.py \
  --preset standard \
  --formats dense,csr,csc \
  --workloads core \
  --repeats 5 \
  --label baseline \
  --output benchmarks/results/baseline.json
```

The `core` workloads are `resolve_x`, `tidy_x`, and `aggregate_x`. Use
`--workloads all` to add layer, raw, mixed-source, embedding, grammar, and
high-level plot-preparation cases. These include `plot_highest_expr_prepare`,
which measures the sparse-aware whole-matrix ranking used by
`plot_highest_expr_genes`:

```bash
PYTHONPATH=src python benchmarks/run_benchmarks.py \
  --preset standard \
  --formats dense,csr,csc \
  --workloads all \
  --repeats 3 \
  --label baseline-full \
  --output benchmarks/results/baseline-full.json
```

Rendering is excluded by default. Add bounded embedding and dotplot renders with:

```bash
PYTHONPATH=src python benchmarks/run_benchmarks.py \
  --preset standard \
  --formats csr \
  --workloads plot_embedding_prepare,plot_dotplot_prepare \
  --include-rendering \
  --repeats 3 \
  --label rendering \
  --output benchmarks/results/rendering.json
```

`extended` uses a 6,000 x 3,000 dense fixture and 20,000 x 10,000 sparse
fixtures. It is opt-in because unoptimized whole-matrix paths can take several
minutes:

```bash
PYTHONPATH=src python benchmarks/run_benchmarks.py \
  --preset extended \
  --formats csr,csc \
  --workloads core \
  --repeats 3 \
  --label extended \
  --output benchmarks/results/extended.json
```

## Scale one dimension

Override `n_obs`, `n_vars`, requested `n_genes`, or `n_groups` with a positive
integer. A comma-separated list creates one case per value, while the remaining
dimensions stay fixed:

```bash
PYTHONPATH=src python benchmarks/run_benchmarks.py \
  --preset smoke \
  --formats csr \
  --workloads resolve_x,tidy_x,aggregate_x \
  --n-obs 250,500,1000 \
  --repeats 3 \
  --label scale-cells \
  --output benchmarks/results/scale-cells.json
```

Only one dimension may contain multiple values in a run, which keeps the sweep
one-factor-at-a-time. Other dimensions may have singleton overrides, for example
`--n-vars 1000 --n-genes 8`. Scaled case IDs include every value that differs
from the preset, such as `smoke[n_obs=250]/csr/resolve_x`; unchanged preset cases
retain their existing IDs.

List individual workload names with:

```bash
python benchmarks/run_benchmarks.py --list-workloads
```

## Compare revisions

Run the same command and seed before and after a change. Then compare the JSON
documents:

```bash
python benchmarks/compare_results.py \
  benchmarks/results/baseline.json \
  benchmarks/results/optimized.json \
  --output benchmarks/results/comparison.md
```

The comparison reports cold time, repeated median time, sampled peak RSS,
retained RSS after garbage collection, and output fingerprints. It can enforce a
regression threshold and unchanged prepared output:

```bash
python benchmarks/compare_results.py \
  benchmarks/results/baseline.json \
  benchmarks/results/optimized.json \
  --fail-regression-pct 5 \
  --fail-on-output-change
```

The comparator rejects different fixtures, seeds, repeat counts, sampling
intervals, platforms, Python builds, or dependency versions. It also requires
the same case set. Use `--allow-incomparable` only to produce an exploratory
table; the resulting report is marked with the mismatches.

## Metrics

Each JSON result records:

- all repeated timings and their median, minimum, and maximum;
- RSS immediately before the call, sampled peak RSS, RSS while the output is
  retained, and RSS after deletion and garbage collection;
- storage used by `X`, layers, raw, embeddings, `obs`, and `var`;
- deep DataFrame or plot-data size, schema, categorical levels, and a content
  fingerprint;
- untimed `stage_sizes` for prepared data and, in core workloads, the projected
  expression matrix and observation frame;
- the resolved `ggann` package path and a SHA-256 digest of its Python source
  tree, plus the effective BLAS and OpenMP thread settings.

The fingerprint and stage sizes are computed after all timed calls. This keeps
pandas hashing and diagnostic matrix selection out of the timings and their RSS
baselines.

RSS is sampled from `/proc/self/statm` every millisecond by default. Change the
interval with `--rss-interval-ms`. On non-Linux systems the fallback is the
process high-water mark, so retained-memory values are not directly comparable
with Linux results. Sampling can miss allocations shorter than the interval.

Use the same machine, dependency environment, thread settings, preset, seed, and
workload list for before-and-after comparisons. The runner pins BLAS and OpenMP
thread counts to one in child processes unless those variables are already set.

## Check downstream artifacts

Snapshot prepared data, mappings, and four rendered figures from
`pbmc68k_reduced`, then compare revisions with numeric tolerances:

```bash
PYTHONPATH=/path/to/baseline python benchmarks/check_invariance.py snapshot \
  --output /tmp/ggann-artifacts-baseline
PYTHONPATH=src python benchmarks/check_invariance.py snapshot \
  --output /tmp/ggann-artifacts-candidate
python benchmarks/check_invariance.py compare \
  /tmp/ggann-artifacts-baseline /tmp/ggann-artifacts-candidate
```
