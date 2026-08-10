# Add a custom depth-aware summary with annplyr

Use ggann helpers for standard plots and annplyr directly when a custom figure
needs an explicitly shaped table or a hard materialization budget.
Grammar plots can instead enforce the same boundary with
`gganndata(..., max_matrix_values=...)`.

The worked question is deliberately practical, and it is the caveat every
analyst meets: library-size normalization does not fully remove the sequencing
depth effect. `depth` is a median split of the measured `n_counts`, so the
summary compares mean CST3 expression between shallowly and deeply sequenced
cells within each PBMC3K cluster, then keeps a standard marker dotplot beside
that custom summary.

`adata.ap.to_df(...)` below projects exactly two genes. Its cumulative budget is
`2 × n_obs`; annplyr rejects an over-budget plan before reading any requested
matrix source. The grouped summary requests one gene and groups by both
`louvain` and `depth`:

```python
summary = adata.ap.summarize(
    x={"mean_CST3": ap.mean(ap.col("CST3"))},
    by=["louvain", "depth"],
    max_matrix_values=adata.n_obs,
)

custom = (
    ggplot(summary, aes("louvain", "mean_CST3", fill="depth"))
    + geom_col(position="dodge")
    + ag.theme_publication()
)
```

The workflow densifies only the single bounded sparse result column before
feeding it to plotnine. Wide annplyr exports intentionally preserve pandas
sparse columns; current plotnine position scales require dense columns.

The depth gap is real and one-directional: in CD14+ monocytes, the largest
myeloid cluster, mean CST3 is 2.88 in deeply sequenced cells against 2.06 in
shallow ones. The executable check asserts on that cluster only, because the
rarest clusters hold as few as one cell per depth bin.

The final `plot_dotplot` call shows that direct annplyr use and ggann helpers can
coexist on the same immutable object. A state fingerprint verifies that neither
workflow modifies the input.

## Executed source

```{literalinclude} ../../examples/vignettes/04_annplyr_interop.py
:language: python
:linenos:
```

Prefer `to_df` for one-row-per-cell tables and `to_tidy` for long
observation-by-feature tables. Always select features explicitly for large
objects; use `max_matrix_values=` when exceeding a known memory boundary must be
an error rather than a performance surprise.
