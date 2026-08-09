# Review a large sparse or backed AnnData

Consider a marker review against a read-only atlas file: the object may have
tens of thousands of genes, while the figure needs four. The safe workflow is
to project those genes before conversion and to keep the source object
immutable.

The important scaling rule is to select columns before converting them to a
plotting table. ggann projects requested genes and observation metadata before
annplyr extraction, so work scales with `n_obs × requested_genes`, not the full
number of variables.

The executable workflow covers:

- CSR and CSC expression matrices;
- aggregated dotplots without whole-matrix densification;
- backed dense, CSR, and CSC `.h5ad` files opened read-only;
- projected `.X` reads with `use_raw=False`;
- mutation fingerprints before and after plotting.

The user-facing call is unchanged across in-memory CSR/CSC objects and a backed
`.h5ad` opened with `backed="r"`:

```python
markers = ["CD3D", "NKG7", "MS4A1", "CST3"]

plot = ag.plot_matrixplot(
    backed,
    markers,
    group_by="cell_type",
    use_raw=False,
)
```

Only the requested expression width is prepared. That keeps the biological
intent visible in the call and avoids the common accidental full-matrix
conversion.

The deterministic matrix is intentionally small so the documentation build is
fast. It exercises storage and ownership contracts, not large-scale performance;
the extended benchmark suite supplies that evidence.

Backed plotting still materializes the selected table required by plotnine; it
does not make rendering itself lazy. Use `downsample=` only when a representative
cell subset is scientifically appropriate. Dotplots and matrixplots aggregate
all cells by default, so the fixture verifies their full-population behavior.

## Executed source

```{literalinclude} ../../examples/vignettes/03_sparse_and_backed.py
:language: python
:linenos:
```

For a hard cumulative read limit in a custom extraction, use annplyr's
`max_matrix_values=` as shown in {doc}`annplyr_interop`.
