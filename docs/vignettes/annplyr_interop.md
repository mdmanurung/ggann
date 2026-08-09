# Add a custom condition-aware summary with annplyr

Use ggann helpers for standard plots and annplyr directly when a custom figure
needs an explicitly shaped table or a hard materialization budget.
Grammar plots can instead enforce the same boundary with
`gganndata(..., max_matrix_values=...)`.

The worked question is deliberately practical: compare mean CD3D expression
between control and stimulated cells within each annotated lineage, then keep a
standard marker dotplot beside that custom summary.

`adata.ap.to_df(...)` below projects exactly two genes. Its cumulative budget is
`2 × n_obs`; annplyr rejects an over-budget plan before reading any requested
matrix source. The grouped summary requests one gene and groups by both
`cell_type` and `condition`:

```python
summary = adata.ap.summarize(
    x={"mean_CD3D": ap.mean(ap.col("CD3D"))},
    by=["cell_type", "condition"],
    max_matrix_values=adata.n_obs,
)

custom = (
    ggplot(summary, aes("cell_type", "mean_CD3D", fill="condition"))
    + geom_col(position="dodge")
    + ag.theme_publication()
)
```

The workflow densifies only the single bounded sparse result column before
feeding it to plotnine. Wide annplyr exports intentionally preserve pandas
sparse columns; current plotnine position scales require dense columns.

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
