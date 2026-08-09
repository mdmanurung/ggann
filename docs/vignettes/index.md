# Executable vignettes

These six workflows run in isolated Python processes during every Sphinx
build. They use the deterministic fixture in `examples/vignettes/_fixture.py`,
set Matplotlib's headless backend, write only to temporary directories, and do
not access the network.

```{toctree}
:maxdepth: 1

scanpy_migration
grammar_of_graphics
sparse_backed
annplyr_interop
publication_panels
scanpy_tradeoffs
```

Run them directly from a source checkout:

```bash
for script in examples/vignettes/[0-9]*.py; do
  python "$script"
done
```

The warning-as-error documentation build runs the same scripts automatically:

```bash
GGANN_DOCS_OFFLINE=1 sphinx-build -W --keep-going -b html docs docs/_build/html
```
