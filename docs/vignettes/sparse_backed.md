# Large sparse and backed AnnData

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

The deterministic matrix is intentionally small so the documentation build is
fast. It exercises storage and ownership contracts, not large-scale performance;
the extended benchmark suite supplies that evidence.

Backed plotting still materializes the selected table required by plotnine. It
does not make rendering itself lazy. Use `downsample=` only when a representative
cell subset is scientifically appropriate; aggregated plots otherwise use all
cells.

## Executed source

```{literalinclude} ../../examples/vignettes/03_sparse_and_backed.py
:language: python
:linenos:
```

For a hard cumulative read limit in a custom extraction, use annplyr's
`max_matrix_values=` as shown in {doc}`annplyr_interop`.
