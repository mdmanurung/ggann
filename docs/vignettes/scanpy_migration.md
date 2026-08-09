# Migrate from Scanpy plotting

Most migrations replace a `sc.pl` call with the corresponding `ggann` helper,
then refine the returned plot with plotnine.

| Scanpy | ggann |
|---|---|
| `sc.pl.embedding(adata, "umap", color=...)` | `ag.plot_embedding(adata, "umap", color=...)` |
| `sc.pl.dotplot(adata, genes, groupby=...)` | `ag.plot_dotplot(adata, genes, group_by=...)` |
| `sc.pl.matrixplot(adata, genes, groupby=...)` | `ag.plot_matrixplot(adata, genes, group_by=...)` |
| `sc.pl.violin(adata, genes, groupby=...)` | `ag.plot_violin(adata, genes, group_by=...)` |

The spelling change from Scanpy's `groupby` to `group_by` is intentional and is
consistent across ggann helpers. `show=` and `return_fig=` are unnecessary:
ggann returns the plot object without drawing it.

```python
scanpy_plot = sc.pl.embedding(
    adata,
    basis="umap",
    color="cell_type",
    show=False,
)

ggann_plot = ag.plot_embedding(
    adata,
    basis="umap",
    color="cell_type",
)
```

Both calls use the same source-selection convention: `.raw` when present,
otherwise `.X`, unless `use_raw=` or `layer=` selects a source explicitly.
Visual parity means the same cells, variables, groups, and statistics—not pixel
identity between Matplotlib and plotnine. See {doc}`../scanpy_parity` for the
full mapping and intentional differences.

## Executed source

```{literalinclude} ../../examples/vignettes/01_scanpy_migration.py
:language: python
:linenos:
```
