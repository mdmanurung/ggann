# Contributing

Bug reports and focused pull requests are welcome. For behavior changes, open
an issue first so the public API and Scanpy-parity implications can be agreed
before implementation.

## Development setup

ggann requires Python 3.12 or newer. Create an isolated environment and install
the development dependencies:

```bash
uv venv .venv --python 3.12
uv pip install --python .venv/bin/python -e ".[dev]"
```

plotnine resolves to 0.16.0a15 or newer within `<0.17` without `--pre`. If you
reuse an old environment with plotnine 0.15 or 0.16.0a14, upgrade it.

annplyr 0.5 and plotnine-extra 0.3.1 are not on PyPI yet. Install them before
the development extras, either from the commits CI pins:

```bash
uv pip install --python .venv/bin/python \
  "annplyr @ git+https://github.com/mdmanurung/annplyr@bb80d1285d41a1e0470311bb887018433c5acc80" \
  "plotnine-extra @ git+https://github.com/mdmanurung/plotnine-extra@b24b997749b67ea147a1f759df3d6bc997f511a1"
uv pip install --python .venv/bin/python -e ".[dev]"
```

or from local sibling checkouts:

```bash
uv pip install --python .venv/bin/python -e ../annplyr -e ../plotnine-extra
uv pip install --python .venv/bin/python -e ".[dev]"
```

The sibling repositories are dependencies, not part of the ggann change set;
do not modify them from a ggann pull request.

## Required checks

The test suite and the documentation build both execute the vignettes, which
read the PBMC3K dataset from a local cache rather than downloading it. Fetch it
once, while online:

```bash
python scripts/fetch_datasets.py
```

Run the same deterministic checks used in continuous integration:

```bash
python scripts/run_pyright.py
ruff check src tests benchmarks docs/extensions examples scripts
ruff format --check src tests benchmarks docs/extensions examples scripts
PYTHONPATH=src python -m pytest -q   # about 7 minutes
GGANN_DOCS_OFFLINE=1 sphinx-build -W --keep-going -b html docs docs/_build/html
python -m build
python -m twine check dist/*
```

### Plotnine coupling

`python -m pytest -m canary -q` runs the plotnine canary tests alone. They assert
effects (compose widths, heights and gap), because some plotnine 0.16.0a15 breaks
were silent. They also run in the normal test job, and plotnine-extra's
`Nightly canary` workflow runs them every night against the newest plotnine
pre-release.

Import private plotnine names only through `src/ggann/_plotnine_compat.py`
(`Environment`, `BUILD_CACHE_ATTRS`, `sub_gridspec`, `grid_shape`,
`figure_edge_fractions`, `pin_grid_spacing`, `missing_coupling()`). After a
plotnine bump, run `pytest -m canary` first; a failure points at that module.
`src/ggann/_compat.py` is a different module, for Scanpy keyword aliases.

Changes to extraction, aggregation, or plotting must also run the benchmark
smoke suite and compare representative prepared data and rendered artifacts.
Include the exact commands and environment versions in the pull request.

## API changes

- Keep helper arguments consistent across plot families.
- Preserve ordinary `plotnine.ggplot` return values for grammar-native plots.
- Add a deprecation period and migration note for unavoidable public changes.
- Add contract tests for signatures, defaults, errors, return types, and input
  ownership.

Do not commit generated build directories, benchmark machine artifacts, or
downloaded datasets.
