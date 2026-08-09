# ggann release readiness

Assessment date: 2026-08-09

## Final verdict: NOT READY

The repository is a locally verified release candidate for documentation, API
contracts, annplyr interoperability, tests, and distribution mechanics. It is
**not ready for publication** because the required annplyr and plotnine-extra
versions are not available from the public package index and the measured
large-sparse Scanpy performance gates fail decisively. No "blazing fast" claim
is supported.

## Status summary

| Area | Result | Evidence |
|---|:---:|---|
| Documentation | **PASS** | Warning-as-error offline Sphinx build; all six deterministic vignettes execute, including the evidence-linked Scanpy trade-off comparison; first-party public docstrings include parameters, returns, errors, notes, and examples. |
| ggplot2-style API | **PASS** | Ordinary `plotnine.ggplot` grammar objects, helper/grammar parity, composability, selectors, signatures, defaults, return boundaries, and representative errors are contract-tested. |
| annplyr compatibility | **PASS** | Public annplyr 0.3 projected extraction; dense/CSR/CSC, views, duplicate/reordered names, X/layer/raw/obs/obsm, backed dense/CSR/CSC, budgets, typed errors, missing values, and no mutation are tested. |
| Performance | **FAIL** | Profiling-guided changes materially improve embedding and matrixplot, but all primary seven-repeat extended-CSR preparation, end-to-end, and peak-RSS gates still fail; the end-to-end geometric speedup is 0.526x. |
| Packaging | **FAIL** | Wheel/sdist mechanics pass, but an index-only install cannot resolve `annplyr>=0.3,<0.4` or `plotnine-extra>=0.3.1,<0.4`. |
| CI configuration | **PASS** | Python 3.12/3.13/3.14, lint/format, tests, docs, smoke benchmark, build/twine, extras import, and clean-wheel jobs are configured; equivalent local commands pass. Hosted workflows were not triggered because no commit was pushed. |

## Provenance

| Item | Value |
|---|---|
| ggann baseline commit | `d397b2ab2d580e8733962aca8ccf77cf56d7e01c` |
| Branch | `agent/performance-api-release-refactor` |
| Candidate version | `0.1.0` |
| Frozen pre-optimization source SHA-256 | `b8c9eb82988546b73b23c3bd5649c1756537eb07c766ebd79a87051227a737b7` |
| Candidate source SHA-256 used by the final Scanpy run | `7e059e89a18746625fb9e3727e2b1f1d7d1ae22050d2b5dfa711832253bf4e0d` |
| Python benchmark runtime | `3.12.13` |
| Platform | Linux 4.18, glibc 2.28, Intel Xeon E5-2690 v3, 24 logical CPUs |
| annplyr baseline tested | `0.3.0`, commit `ef136845cc425eb157530247a626eac4d23a2219` |
| annplyr v0.3 implementation snapshot tested | `0.3.0`, commit `a6edfb2c05be3f865fa46d77b2be22db454fec2d` |
| annplyr current sibling tested separately | `0.3.0`, commit `cadb0e77eeb4ce9eef6a47deb9a51b93e6abf617` |
| plotnine-extra snapshot tested | `0.3.1`, commit `fcbd88f64760a2dd6e8200028a1a00909aa2c9f8` |

The ggann worktree was clean at the start. `../annplyr` was already dirty and
changed commits externally during this work; ggann did not modify it. Clean
immutable archives of the annplyr baseline, implementation, and current release
commits were used in isolated environments, and the final current-sibling
compatibility test was run separately. The current `cadb0e7` commit changes only
release documentation/configuration relative to `a6edfb2`; package code is
identical. The local
`../plotnine-extra` checkout remained untouched.

Benchmark package versions were Scanpy 1.12.3, AnnData 0.13.2, annplyr 0.3.0,
NumPy 2.4.6, pandas 3.0.5, SciPy 1.18.0, plotnine 0.15.7,
plotnine-extra 0.3.1, and Matplotlib 3.11.1. BLAS/OpenMP thread counts were one;
the renderer was Agg at 6 x 4.5 inches, 80 DPI, PNG.

## What changed

- Runtime Git dependencies were replaced by bounded version ranges.
- AnnData-to-table reads now delegate to public annplyr 0.3 projected APIs.
- Mixed X/raw/layer/obsm grammar requests receive a cumulative pre-read
  `max_matrix_values` check; annplyr typed errors propagate.
- Duplicate/reordered observation names are assembled positionally, and
  pandas-backed `obsm` columns resolve by their actual positional labels.
- Group means and detection fractions now use one projected annplyr read and a
  sparse-native vectorized reduction. This removed the repeated
  `annplyr.summarize` hotspot.
- Compatible expression/grouping and embedding/colour/facet requests are now
  fused into one public annplyr extraction call, while cumulative budgets,
  typed errors, and positional reconstruction remain intact.
- Compact group-by-gene frames are built without pandas `melt`, fresh plotnine
  objects/themes avoid repeated internal deep copies, and embedding points skip
  a redundant same-colour outline while preserving their displayed diameter.
- Sparse expression `NaN` values now match pandas/Scanpy semantics: means skip
  them and detection fractions treat them as not detected.
- README, concepts, quickstart, parity/migration/performance guides, five
  core executable vignettes, full public docstrings, release metadata, and CI were
  added or hardened.
- A sixth executable vignette now renders four matched ggann/Scanpy pairs,
  validates their prepared payloads, demonstrates grammar advantages, and loads
  the committed benchmark JSON rather than running the expensive suite in CI.

## Public API changes

All changes are additive or documented consistency fixes; no public symbol was
removed and return types are unchanged.

The profiling follow-up introduced no new public signatures or return-type
changes. Existing calls require no migration; the table below remains the
complete release-candidate API delta.

| Before | After | Migration |
|---|---|---|
| `gganndata(adata, mapping, *, layer=None, use_raw=None, add_theme=True)` | Adds keyword-only `max_matrix_values=None` | Existing calls are unchanged. Set an integer to reject a mapping before extraction when cumulative expression/obsm values exceed the budget. |
| `embedding_coords(adata, basis, n=2)` | Adds keyword-only `max_matrix_values=None` | Existing calls are unchanged. The budget is `n_obs * selected_coordinates`. |
| Grouped summaries could retain a row whose grouping key was missing | Rows missing any `group_by`/`split_by` key are omitted | Fill missing labels explicitly before plotting when they should be a category. |
| Sparse grouped expression-`NaN` behavior depended on the backend reduction | Means skip expression `NaN`; fractions count it as not detected | Finite matrices are unchanged. |
| Float32 grouped accumulation | Float64 grouped accumulation | Results remain equivalent within `rtol=1e-6`, `atol=1e-7`; no call-site change. |

All plotnine-native helpers still return ordinary `plotnine.ggplot` objects.
`plot_clustermap` and `plot_upset` remain the documented grid-backend
exceptions. The canonical shared names remain `color`, `group_by`, `split_by`,
`layer`, `use_raw`, `downsample`, `random_state`, and `categories_order`.

## annplyr compatibility evidence

| Matrix | Python | AnnData | annplyr | Scanpy | NumPy / pandas / SciPy | plotnine | Result |
|---|---|---|---|---|---|---|---|
| Lowest declared stack | 3.12.13 | 0.12.0 | 0.3.0 (`ef136845`) | 1.11.0 | 1.26.4 / 2.2.0 / 1.12.0 | 0.15.3 | **86 passed** focused extraction/source/API tests |
| Latest compatible | 3.12.13 | 0.13.2 | 0.3.0 (`cadb0e7`) | 1.12.3 | 2.4.6 / 3.0.5 / 1.18.0 | 0.15.7 | **295 passed** |
| Latest compatible | 3.13.11 | 0.13.2 | 0.3.0 (`ef136845`) | 1.12.3 | 2.4.6 / 3.0.5 / 1.18.0 | 0.15.7 | **295 passed** |
| Latest compatible | 3.14.5 | 0.13.2 | 0.3.0 (`ef136845`) | 1.12.3 | 2.4.6 / 3.0.5 / 1.18.0 | 0.15.7 | **295 passed** |
| Current local sibling | 3.12.13 | 0.13.2 | 0.3.0 (`cadb0e7`) | 1.12.3 | 2.4.6 / 3.0.5 / 1.18.0 | 0.15.7 | **86 passed** focused compatibility tests |

The `annplyr>=0.3,<0.4` requirement accurately implies
`anndata>=0.12,<1`. annplyr has no public index release, so “lowest”, “stable”,
and “newest” cannot yet be represented as distinct published versions. The
immutable v0.3 baseline and current-release commits above form the provisional
compatibility matrix.

Covered contracts include:

- dense NumPy, pandas-backed matrices, CSR/CSC matrices and arrays;
- AnnData views and independent objects;
- `.X`, named layers, `.raw`, `obs`, `var`, and `obsm` resolution;
- backed dense/CSR/CSC projections and reordered backed views;
- category ordering, missing values, duplicate names, reordered rows, sparse
  preservation, and expression `NaN` values;
- projection before conversion, cumulative grammar budgets, exact-budget
  success, over-budget zero-read rejection, and typed `AnnplyrError` propagation;
- input ownership and no accidental AnnData mutation.

High-level plotting helpers project requested width but do not each expose a
hard `max_matrix_values` argument. Use the grammar boundary or direct annplyr
extraction when a hard budget is required; this limitation is explicit in the
concepts and large-data guides.

## Matched Scanpy benchmark

The primary evidence is versioned at
`benchmarks/results/scanpy-extended-csr.json`. It uses one immutable 20,000 x
10,000 CSR AnnData, 32 requested genes, 32 groups, seven repeated warm samples,
separate cold samples, fresh-process end-to-end RSS, and numerical payload
validation before timings are accepted. Times are median with `[min, max]`;
speedup is Scanpy time divided by ggann time. Construction is diagnostic because
the public APIs defer different amounts of artist work.

| Workload | Prep: ggann vs Scanpy | Prep speedup | Construction: ggann vs Scanpy | Render: ggann vs Scanpy | End to end: ggann vs Scanpy | E2E speedup | Peak / retained RSS: ggann vs Scanpy |
|---|---:|---:|---:|---:|---:|---:|---:|
| Embedding, categorical colour | 13.3 `[13.0, 14.2]` vs 4.3 `[4.2, 4.4]` ms | 0.33x | 14.6 vs 94.4 ms | 382.6 vs 313.4 ms | 770.4 `[763.6, 796.3]` vs 426.7 `[422.2, 442.1]` ms | 0.55x | 29.09 / 28.17 vs 5.60 / 5.60 MiB |
| Dotplot | 49.5 `[48.3, 59.6]` vs 17.2 `[16.1, 22.4]` ms | 0.35x | 53.5 vs 22.9 ms | 457.8 vs 213.0 ms | 738.0 `[727.7, 768.0]` vs 367.9 `[362.2, 387.6]` ms | 0.50x | 17.00 / 17.00 vs 12.12 / 12.12 MiB |
| Matrixplot | 51.6 `[48.9, 67.4]` vs 12.0 `[11.3, 14.9]` ms | 0.23x | 55.7 vs 15.2 ms | 413.9 vs 217.1 ms | 656.7 `[648.2, 671.6]` vs 346.5 `[336.8, 355.7]` ms | 0.53x | 17.16 / 17.16 vs 7.73 / 7.73 MiB |

| Gate | Required | Measured | Result |
|---|---:|---:|:---:|
| Embedding preparation | >=2.00x | 0.33x | **FAIL** |
| Dotplot preparation | >=2.00x | 0.35x | **FAIL** |
| Matrixplot preparation | >=2.00x | 0.23x | **FAIL** |
| Every primary end to end | >=0.91x | 0.50x to 0.55x | **FAIL** |
| End-to-end geometric mean | >1.00x | 0.526x | **FAIL** |
| Large-sparse peak RSS | <=1.10x Scanpy | 1.40x to 5.19x | **FAIL** |

The full repeated dense/CSR/CSC by `.X`/layer/raw large matrix and a real-PBMC
cross-library benchmark are not complete. This incompleteness is an additional
performance blocker, not a reason to soften the measured failures.

## Frozen ggann baseline comparison

The exact pre-optimization candidate was rerun from its frozen source snapshot
with the same environment, fixture, seed, and seven repetitions. Raw documents
are versioned as `scanpy-extended-csr-before.json` and
`scanpy-extended-csr.json`.

| Workload | Prep before → after | Construction before → after | Render before → after | End to end before → after | E2E peak RSS before → after |
|---|---:|---:|---:|---:|---:|
| Embedding | 16.6 → 13.3 ms (-19.9%) | 23.1 → 14.6 ms (-36.7%) | 540.5 → 382.6 ms (-29.2%) | 940.2 → 770.4 ms (-18.1%) | 28.70 → 29.09 MiB (+1.4%) |
| Dotplot | 50.7 → 49.5 ms (-2.5%) | 61.7 → 53.5 ms (-13.2%) | 465.0 → 457.8 ms (-1.5%) | 745.6 → 738.0 ms (-1.0%) | 16.90 → 17.00 MiB (+0.6%) |
| Matrixplot | 59.1 → 51.6 ms (-12.7%) | 69.8 → 55.7 ms (-20.2%) | 401.6 → 413.9 ms (+3.0%) | 666.6 → 656.7 ms (-1.5%) | 16.19 → 17.16 MiB (+6.0%) |

Every primary preparation and end-to-end timing passes the 5% self-regression
limit. The single fresh-child matrixplot RSS increase exceeds it; because one
probe cannot distinguish allocator noise from a real ~1 MiB retention change,
it remains an explicit blocker rather than being accepted without evidence.

The concrete diagnostic profiler independently records:

| Profiled path | Before | After | Change |
|---|---:|---:|---:|
| Embedding preparation | 20.3 ms | 14.5 ms | -28.5% |
| Embedding end to end | 961.2 ms | 777.1 ms | -19.1% |
| Dotplot preparation | 56.2 ms | 52.3 ms | -6.8% |
| Dotplot end to end | 758.1 ms | 734.3 ms | -3.1% |
| Matrixplot preparation | 65.4 ms | 53.5 ms | -18.1% |
| Matrixplot end to end | 685.1 ms | 636.6 ms | -7.1% |

The profiler JSON and raw timing samples are versioned in
`benchmarks/results/profile-hotspots-before.json` and
`profile-hotspots-after.json`; raw pstats are retained under
`/tmp/ggann-profile-before-hotspots` and `/tmp/ggann-profile-after-hotspots`.

## Downstream invariance

The candidate was compared with the exact baseline on the cached real
`scanpy.datasets.pbmc68k_reduced` dataset:

| Artifact | Result |
|---|:---:|
| Resolved and aggregated data | **PASS**, `rtol=1e-6`, `atol=1e-7` |
| Embedding/dotplot/matrixplot/heatmap prepared data | **PASS** |
| Plot mappings, labels, and scales | **PASS**, exact |
| Dotplot/matrixplot/heatmap PNGs | **PASS**, pixel-identical (`max=0`, `mean=0`) |
| Embedding PNG | **PASS** under the explicit outline-equivalence tolerance: mean normalized error 0.00194 ≤ 0.002 and PSNR 35.27 dB ≥ 35 dB; prepared data, mapping, labels, scales, population, and marker diameter are unchanged |

## Packaging and publication order

Release metadata contains no URL or Git runtime requirements. The final local
distribution verification is recorded below; the packaging category remains
failed because an index-only installation cannot resolve the declared runtime
dependencies.

Required publication order:

1. Publish annplyr 0.3.0 from the intended `cadb0e7` release source.
2. Publish plotnine-extra 0.3.1.
3. Remove the temporary pinned CI/Read the Docs bootstrap checkouts.
4. Rebuild ggann and repeat an index-only clean-wheel install plus `pip check`.
5. Rerun the complete repeated performance matrix and proceed only if every
   gate passes.

## Verification commands and outcomes

| Check | Command | Outcome |
|---|---|---|
| Python 3.12 full suite | `PYTHONPATH=src NUMBA_CACHE_DIR=/tmp/ggann-numba-cache-py312-final MPLCONFIGDIR=/tmp/ggann-mpl-cache-py312-final /tmp/ggann-wheel-final-KtBMp5/venv/bin/python -m pytest -q` | **295 passed**, 16 expected warnings |
| Python 3.13 full suite | `PYTHONPATH=src NUMBA_CACHE_DIR=/tmp/ggann-numba-cache-py313 MPLCONFIGDIR=/tmp/ggann-mpl-cache-py313 /tmp/ggann-matrix-py313-2U4sZk/bin/python -m pytest -q` | **295 passed**, 16 expected warnings |
| Python 3.14 full suite | `PYTHONPATH=src NUMBA_CACHE_DIR=/tmp/ggann-numba-cache-py314 MPLCONFIGDIR=/tmp/ggann-mpl-cache-py314 /tmp/ggann-matrix-py314-pHFqlv/bin/python -m pytest -q` | **295 passed**, 16 expected warnings |
| Lowest dependency stack | `PYTHONPATH=src /tmp/ggann-matrix-oldest-py312-LWFy4e/bin/python -m pytest -q tests/test_expression_projection.py tests/test_resolve.py tests/test_layers.py tests/test_prefixes.py tests/test_api_consistency.py tests/test_release_metadata.py tests/test_aggregate.py` | **86 passed** |
| Current sibling annplyr | Same focused command after installing an immutable archive of read-only `../annplyr` at `cadb0e7` | **86 passed** |
| Documentation | `PYTHONPATH=src GGANN_DOCS_OFFLINE=1 MPLCONFIGDIR=/tmp/ggann-mpl-docs-profile /tmp/ggann-wheel-final-KtBMp5/venv/bin/sphinx-build -W --keep-going -b html docs /tmp/ggann-docs-after-profile` | **PASS**; all six vignettes executed offline |
| Lint and format | `ruff check src tests benchmarks docs/extensions examples` and `ruff format --check src tests benchmarks docs/extensions examples` | **PASS** |
| YAML/CFF syntax | PyYAML `BaseLoader` over both workflows, `.readthedocs.yaml`, and `CITATION.cff` | **PASS** |
| Hotspot profiler smoke | `PYTHONPATH=src /tmp/ggann-wheel-final-KtBMp5/venv/bin/python benchmarks/profile_hotspots.py --preset smoke --format csr --workloads embedding_categorical,dotplot,matrixplot --stages preparation --repeats 1 --output /tmp/ggann-profile-final-smoke.json --profile-dir /tmp/ggann-profile-final-smoke` | **PASS**; stage timings, allocations, RSS, frame sizes, sparse conversions, extraction calls, and pstats retained |
| Real-data invariance | `python benchmarks/check_invariance.py compare /tmp/ggann-artifacts-candidate-final-formatted /tmp/ggann-artifacts-after-profile-opt --allow-image-difference embedding --image-mean-tolerance 0.002 --image-psnr-min 35 --report /tmp/ggann-profile-opt-invariance-tolerated.md` | **PASS**; prepared values/mappings/labels/scales exact or within stated numeric tolerance, three PNGs pixel-identical, embedding within explicit antialiasing tolerance |
| Repeated Scanpy comparison | `PYTHONPATH=src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python benchmarks/compare_scanpy.py --preset extended --formats csr --workloads primary --sources x --repeats 7 --seed 20260809 --include-cold-start --isolated-memory-stages preparation,end_to_end --output /tmp/ggann-final-extended-7.json --report /tmp/ggann-final-extended-7.md` | Completed with payload equivalence and input immutability; JSON gate status **FAIL** |
| Frozen baseline comparison | Repeat the command above with `PYTHONPATH=/tmp/ggann-profile-before-eICs5m/src`, writing `/tmp/ggann-final-baseline-extended-7.json`; compare it with the candidate's identically configured raw samples | **PASS** for the 5% preparation/end-to-end timing gate; matrixplot fresh-child peak RSS increased 6.0% and remains a blocker |
| Wheel and sdist | `/tmp/ggann-release-py312/bin/python -m build --outdir /tmp/ggann-dist-20260809-push` followed by `/tmp/ggann-release-py312/bin/python -m twine check /tmp/ggann-dist-20260809-push/*` | **PASS**; final wheel and sdist built from this report revision; required benchmark/vignette/report data included in the sdist and excluded appropriately from the wheel |
| Clean wheel | Install the built wheel plus immutable annplyr 0.3.0 and plotnine-extra `fcbd88f` snapshots into `/tmp/ggann-cleanwheel-profile`, then import all extras, render, exercise the budget boundary, and run `pip check` from `/tmp` | **PASS**; import came from site-packages, 7,902-byte PNG, ordinary `plotnine.ggplot`, typed budget error, all extras imported, no broken requirements |
| Public-index-only install | `uv pip install <wheel>` without dependency snapshots | **FAIL** as expected: required annplyr and plotnine-extra versions unavailable |

## Remaining blockers and unsupported cases

- Every measured primary Scanpy speed and memory gate fails.
- Plotnine scale/guide/theme/axis work and Matplotlib artist construction remain
  the dominant render costs; annplyr's projected-frame/name-position adapters
  remain the dominant preparation cost after extraction-call fusion.
- The full repeated large dense/CSR/CSC and X/layer/raw matrix is incomplete.
- `plot_highest_expr_genes` has an unresolved zero-total-observation semantic
  difference in the cross-library comparator.
- Public-index installation is impossible until annplyr 0.3.0 and
  plotnine-extra 0.3.1 are published.
- Matrixplot fresh-child end-to-end peak RSS is 6.0% above the frozen ggann
  candidate in the single repeated comparison, despite invariant outputs.
- Hosted CI has not run because this task intentionally did not commit or push.

These blockers require the verdict **NOT READY**.
