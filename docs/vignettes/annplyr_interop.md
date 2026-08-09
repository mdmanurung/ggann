# annplyr interoperability

Use ggann helpers for standard plots and annplyr directly when a custom figure
needs an explicitly shaped table or a hard materialization budget.
Grammar plots can instead enforce the same boundary with
`gganndata(..., max_matrix_values=...)`.

`adata.ap.to_df(...)` below projects exactly two genes. Its cumulative budget is
`2 × n_obs`; annplyr rejects an over-budget plan before reading any requested
matrix source. The grouped summary uses the same projection machinery, then
densifies its single bounded sparse result column before feeding the table to
plotnine with `theme_ggann()`. Wide annplyr exports intentionally preserve
pandas sparse columns; current plotnine position scales require dense columns.

The final helper call shows that direct annplyr use and ggann use can coexist on
the same immutable object. A state fingerprint verifies that neither workflow
modifies the input.

## Executed source

```{literalinclude} ../../examples/vignettes/04_annplyr_interop.py
:language: python
:linenos:
```

Prefer `to_df` for one-row-per-cell tables and `to_tidy` for long
observation-by-feature tables. Always select features explicitly for large
objects; use `max_matrix_values=` when exceeding a known memory boundary must be
an error rather than a performance surprise.
