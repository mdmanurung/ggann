# ggann 0.1.0 verification record

Assessment date: 2026-11-02. Branch `release/0.1.0-readiness`.

The package source is ready for release. Publication to PyPI is blocked on
two upstream packages, `annplyr` and `plotnine-extra`, which must reach PyPI
first (see [Release blockers](#release-blockers)).

## Outcome

| Area | Result | Evidence |
|---|:---:|---|
| Test suite | **PASS** | 680 tests, Python 3.12.13, `pytest -n 8` |
| Backlog | **PASS** | All 96 evaluation rows completed, documented or deferred with a reason |
| Static quality | **PASS** | ruff check and format clean; pyright 0 errors over `src/ggann` |
| Documentation | **PASS** | Offline `sphinx-build -W --keep-going` exits 0 and executes every vignette |
| Scanpy agreement | **PASS** | Dotplot/matrixplot means and fractions match Scanpy within `atol=rtol=1e-6` (max 3e-7) |
| Statistics | **PASS** | `stat_cor`, `stat_regline_equation`, `stat_compare_means` match SciPy on the same inputs |
| Distribution | **PASS** | sdist and wheel build; `twine check --strict`; wheel imports outside the checkout; `pip check` clean |
| Public install | **BLOCKED** | `annplyr` 404 and `plotnine-extra` only 0.2.0 on PyPI |

## Runtime

Python 3.12.13 in `.venv-rc`, locked in `plans/evaluation/env-rc.lock`.
First-party dependencies come from pinned commits, as in CI: annplyr
`cadb0e77` (0.3.0) and plotnine-extra `fcbd88f6` (0.3.1).

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

Partial and deferred rows are listed for users in `docs/limitations.md`.
They include integer-data sina jitter (upstream plotnine), clustermap
clipping at publication sizes, no helper-output late grammar fields, no
`dendrogram=` on dotplot/matrixplot, negative `n_genes`, and fonts verified
only against a simulated Arial. Seven upstream drafts are in
`plans/evaluation/upstream/`; none were filed.

## Performance

Numbers measured on `res-hpc-exe048` (Xeon E5-2690 v3, 8 of 24 CPUs) on
2026-11-01, with other jobs on the node. Details are in `docs/performance.md`.

- Tracksplot: 700 cells 6.12 s to 0.52 s; 2,800 cells 22.89 s to 0.99 s.
- Embedding point density: exact KDE stays the default and warns above
  20,000 cells. Opt-in `method="binned"` takes 0.008 s at 20,000 cells against
  27 s exact, with Spearman 1.0000 and max normalized error 0.016.
- The exact-method target (under 5 s at 30,000 cells) is still unmet.

## Release blockers

`pip install ggann` cannot resolve until both first-party dependencies are on
PyPI at the declared versions. In order:

1. Publish `annplyr` 0.3.x to PyPI.
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
performance-regression gate on a 20,000 × 10,000 fixture. Those checks predate
the changes above and were not re-run as such; deliberate default changes
(category order, palettes, stacked violin) now alter several of the compared
images. The raw results stay in `benchmarks/results/`.
