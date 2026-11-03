# ggann 0.1.0 verification record

Assessment date: 2026-11-02 (gate re-run included). Branch `release/0.1.0-readiness`.

The package source is ready for release. Publication to PyPI is blocked on
two upstream packages, `annplyr` and `plotnine-extra`, which must reach PyPI
first (see [Release blockers](#release-blockers)).

## Outcome

| Area | Result | Evidence |
|---|:---:|---|
| Test suite | **PASS** | 755 tests, Python 3.12.13, `python -m pytest -q` |
| Backlog | **PASS** | All 96 evaluation rows completed, documented or deferred with a reason |
| Static quality | **PASS** | ruff check and format clean; pyright 0 errors over `src/ggann` |
| Documentation | **PASS** | Offline `sphinx-build -W --keep-going` exits 0 and executes every vignette |
| Scanpy agreement | **PASS** | Dotplot/matrixplot means and fractions match Scanpy within `atol=rtol=1e-6` (max 3e-7) |
| Statistics | **PASS** | `stat_cor`, `stat_regline_equation`, `stat_compare_means` match SciPy on the same inputs |
| Distribution | **PASS** | sdist and wheel build; `twine check --strict`; wheel imports outside the checkout; `pip check` clean |
| Self-regression | **PASS** (claim stages) | Preparation and end-to-end within 5% of frozen `30e3359` across two counterbalanced passes; embedding construction +12% (diagnostic stage), explained below |
| Public install | **BLOCKED** | `annplyr` 404 and `plotnine-extra` only 0.2.0 on PyPI |

## Runtime

Python 3.12.13 in `.venv-rc`, locked in `plans/evaluation/env-rc.lock`.
First-party dependencies come from pinned commits, as in CI: annplyr
`ae2db54` (0.5.0) and plotnine-extra `fcbd88f6` (0.3.1).

| Package | Version |
|---|---|
| anndata | 0.13.4 |
| scanpy | 1.12.5 |
| plotnine | 0.15.8 |
| pandas | 3.1.0 |
| numpy | 2.5.4 |
| scipy | 1.18.1 |
| matplotlib | 3.11.2 |
| statsmodels | 0.15.0 |
| decoupler | 2.3.1 |
| PyComplexHeatmap | 1.8.5 |
| marsilea | 0.8.1 |
| pynebulosa | 0.1.1 |

CI still runs Python 3.12, 3.13 and 3.14. This record covers 3.12 only; the
other two run in CI on the next push.

## What changed since the 2026-08-09 record

A six-report evaluation found 96 backlog rows. The work is recorded per area
in `plans/evaluation/changes/` and routed from
`plans/evaluation/remaining_todos.md`. Binding decisions are in
`plans/evaluation/decisions.md`. `CHANGELOG.md` lists every user-visible
change, with figure-altering changes first.

Changes that affect scientific results:

- Cells with a missing `group_by` or `split_by` are dropped from every grouped
  plot and statistic, with a warning that gives the count. Proportions use the
  kept cells as the denominator. Before, per-cell plots kept them and some
  helpers crashed.
- Helper requests are validated at the call. A bare gene string is one gene;
  unknown fields raise instead of being skipped.
- `plot_ma` and `plot_volcano` share one significance rule:
  `padj < cutoff` and `|LFC| >= lfc`.
- `pseudobulk(use_raw=True)` reads counts and genes from `adata.raw` when the
  gene sets differ.
- Rank-genes `values_to_plot` draws the stored DE statistic, never an
  expression mean.
- Jitter and sina spread are reproducible under `random_state`, across
  repeated draws and exports.

Regression tests: one per repaired row, in `tests/test_release_core.py` (15),
`test_release_deqc.py` (42), `test_release_dist.py` (39),
`test_release_emb.py` (38) and `test_release_pub.py` (21). Where a row fixed
a bug, its test was checked to fail on the pre-change source.

## Known limitations

Partial and deferred rows are listed for users in `docs/notes.md`.
They include integer-data sina jitter (upstream plotnine), clustermap
clipping at publication sizes, no helper-output late grammar fields, no
`dendrogram=` on dotplot/matrixplot, negative `n_genes`, and fonts verified
only against a simulated Arial. Seven upstream drafts are in
`plans/evaluation/upstream/`; none were filed.

## Performance

Numbers measured on `res-hpc-exe048` (Xeon E5-2690 v3, 8 of 24 CPUs) on
2026-11-01, with other jobs on the node. Details are in `docs/notes.md`.

- Tracksplot: 700 cells 6.12 s to 0.52 s; 2,800 cells 22.89 s to 0.99 s.
- Embedding point density: exact KDE stays the default and warns above
  20,000 cells. Opt-in `method="binned"` takes 0.008 s at 20,000 cells against
  27 s exact, with Spearman 1.0000 and max normalized error 0.016.
- The exact-method target (under 5 s at 30,000 cells) is still unmet.

### Self-regression gate against frozen `30e3359`

Re-run on 2026-11-02 with the 2026-08-09 settings (synthetic extended CSR
fixture, seed 20260809, 7 repeats, one BLAS thread), both trees in
`.venv-rc`. The node carried load 4 to 13 from other jobs, so the run was
counterbalanced: frozen, candidate, frozen, candidate. Ratios below are the sum
of the two candidate medians over the sum of the two frozen medians.

| Workload | Preparation | Construction | Render | End to end |
|---|---:|---:|---:|---:|
| Embedding | 1.045 | **1.121** | 1.031 | 1.041 |
| Dotplot | 0.996 | 1.035 | 1.016 | 1.015 |
| Matrixplot | 1.024 | 1.036 | 1.006 | 1.022 |

Memory (end-to-end peak and retained RSS) is within 5% in every workload.

The first pass alone failed 7 checks while load rose from 10 to 13; the
reversed pass shows that most of it was drift. Embedding construction is a real
+1.5 ms per plot (+12%) in both passes. A profile traces it to the new
colour-scale path: stored-palette reuse, Scanpy's default palette as the
fallback and numeric-to-continuous dispatch (rows D3 and B5), about 2,200
extra small calls per plot. The harness treats construction as a diagnostic
stage; only preparation and end-to-end support speed claims, and both pass.
Raw results: `benchmarks/results/rerun-2026-11-02-*`.

## Release blockers

`pip install ggann` cannot resolve until both first-party dependencies are on
PyPI at the declared versions. In order:

1. Publish `annplyr` 0.5.0 to PyPI. ggann now requires `annplyr>=0.5,<0.6`
   (the 0.3 pin was widened on 2026-11-02 after porting to annplyr's 0.5 API).
2. Publish `plotnine-extra` 0.3.1 from a commit that exports `DimPlot`,
   `FeatureDimPlot` and `ggvolcano` (not the stale `v0.3.1` tag).
3. Install `ggann[heatmap,density,upset,pseudobulk]` from TestPyPI in a clean
   environment and run one vignette.
4. Remove the `vendor/` checkout steps from `.github/workflows/ci.yml` and
   `docs.yml`; confirm CI on 3.12, 3.13 and 3.14.
5. Set the changelog date, add `date-released` and the Zenodo DOI to
   `CITATION.cff`, tag `v0.1.0`, upload.

## Reproduction commands

```bash
python scripts/fetch_datasets.py
python scripts/run_pyright.py
ruff check src tests benchmarks docs/extensions examples scripts
ruff format --check src tests benchmarks docs/extensions examples scripts
python -m pytest -q -n 8
GGANN_DOCS_OFFLINE=1 MPLBACKEND=Agg \
  sphinx-build -W --keep-going -b html docs docs/_build/html
python -m build --outdir dist
twine check --strict dist/*
```

## Historical record (2026-08-09)

The earlier record reported 386 tests on Python 3.12, 3.13 and 3.14,
pixel-identical legacy rendering against frozen commit `30e3359`, and a 5%
performance-regression gate on a 20,000 × 10,000 fixture. The timing gate was re-run
on 2026-11-02 (see Performance). The pixel comparisons were not: deliberate default changes
(category order, palettes, stacked violin) now alter several of the compared
images. The raw results stay in `benchmarks/results/`.
