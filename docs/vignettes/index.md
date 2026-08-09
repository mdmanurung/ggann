# Executable vignettes

These six workflows run in isolated Python processes during every Sphinx build.
The first-use and publication workflows use Scanpy's bundled
`pbmc68k_reduced` data; the focused grammar, storage, interoperability, and
comparison workflows use a deterministic PBMC-like control/stimulation
fixture. They set Matplotlib's headless backend, write only to temporary
directories by default, and do not access the network.

| Workflow | Practical question | What it demonstrates |
|---|---|---|
| {doc}`scanpy_migration` | How do I turn a familiar Scanpy marker review into reusable figure objects? | Real PBMC UMAP, dotplot, violin, unequal composition, and input immutability |
| {doc}`grammar_of_graphics` | Does a condition-associated MKI67 pattern persist across the embedding? | Explicit AnnData selectors, bounded extraction, facets, scales, and ordinary plotnine layers |
| {doc}`sparse_backed` | Can I review marker programs without loading an atlas-sized file eagerly? | CSR/CSC projection, read-only backed files, and mutation fingerprints |
| {doc}`annplyr_interop` | How do I add a custom condition-aware summary beside standard helpers? | Direct annplyr aggregation, a hard materialization budget, and plotnine composition |
| {doc}`publication_panels` | How do I turn a biological claim into a submission-ready multi-panel figure? | Real PBMC data, shared vocabulary, exact sizing, editable export, and an evidence manifest |
| {doc}`scanpy_tradeoffs` | Where does ggann add ergonomics, and what does it cost? | Matched payloads, rendered outputs, grammar advantages, timing, and memory evidence |

The generated fixture is deliberately small enough for documentation CI, but
its cell types, markers, sparse matrices, conditions, layers, and backed-file
paths represent the decisions made in routine single-cell work. The two real
PBMC workflows show the same API on a bundled biological dataset.

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
