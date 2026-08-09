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

Until annplyr 0.3 and plotnine-extra 0.3.1 are available from the configured
package index, install their local sibling checkouts first:

```bash
uv pip install --python .venv/bin/python -e ../annplyr -e ../plotnine-extra
uv pip install --python .venv/bin/python -e ".[dev]"
```

The sibling repositories are dependencies, not part of the ggann change set;
do not modify them from a ggann pull request.

## Required checks

Run the same deterministic checks used in continuous integration:

```bash
ruff check src tests benchmarks docs/extensions examples
ruff format --check src tests benchmarks docs/extensions examples
pytest -q
GGANN_DOCS_OFFLINE=1 sphinx-build -W --keep-going -b html docs docs/_build/html
python -m build
python -m twine check dist/*
```

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
