# Executable vignettes

These six workflows run in isolated Python processes during every Sphinx build.
All six use the same real dataset: `scanpy.datasets.pbmc3k_processed`, the
2,638-cell 10x Genomics PBMC sample behind Scanpy's and Seurat's clustering
tutorials. Its Louvain cell types and marker genes are the ones a reader who has
followed either tutorial already recognises.

Fetch the dataset once with `python scripts/fetch_datasets.py`. The workflows
then read that cached copy, set Matplotlib's headless backend, write only to
temporary directories by default, and reach no network themselves; with
`GGANN_DOCS_OFFLINE=1` a missing cache is an error rather than a download.

| Workflow | Practical question | What it demonstrates |
|---|---|---|
| {doc}`scanpy_migration` | How do I turn a familiar Scanpy marker review into reusable figure objects? | PBMC3K UMAP, dotplot, violin, unequal composition, and input immutability |
| {doc}`grammar_of_graphics` | Is the NKG7 cytotoxic program confined to the lymphoid compartment? | Explicit AnnData selectors, bounded extraction, facets, scales, and ordinary plotnine layers |
| {doc}`sparse_backed` | Can I review marker programs without loading an atlas-sized file eagerly? | CSR/CSC projection, read-only backed files, and mutation fingerprints |
| {doc}`annplyr_interop` | Does apparent marker expression still track sequencing depth after normalization? | Direct annplyr aggregation, a hard materialization budget, and plotnine composition |
| {doc}`publication_panels` | How do I turn a biological claim into a submission-ready multi-panel figure? | Shared vocabulary, exact sizing, editable export, and an evidence manifest |
| {doc}`scanpy_tradeoffs` | Where does ggann add ergonomics, and what does it cost? | Matched payloads, rendered outputs, grammar advantages, timing, and memory evidence |

PBMC3K is small enough for documentation CI at 2,638 cells and 1,838
highly variable genes, while its cell types, markers, sparse matrices, layers,
and backed-file paths represent the decisions made in routine single-cell work.
Two derived columns appear alongside the published ones: `compartment` recodes
the Louvain labels into lymphoid and myeloid lineages, and `depth` is a median
split of the measured `n_counts`. Both are documented in
`examples/vignettes/_fixture.py`.

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
python scripts/fetch_datasets.py
for script in examples/vignettes/[0-9]*.py; do
  python "$script"
done
```

The warning-as-error documentation build runs the same scripts automatically:

```bash
python scripts/fetch_datasets.py
GGANN_DOCS_OFFLINE=1 sphinx-build -W --keep-going -b html docs docs/_build/html
```
