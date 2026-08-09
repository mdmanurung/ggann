# ggann + annplyr versus Scanpy: choosing the right plotting path

Scanpy and ggann solve different plotting problems. Scanpy's canonical plots
use direct Matplotlib implementations and are currently faster and leaner.
ggann pays for a plotnine grammar so a returned plot can be extended with
ordinary layers, scales, facets, annotations, and themes. annplyr supplies the
explicit, projected AnnData-to-table boundary beneath that grammar.

The examples below execute offline in
`examples/vignettes/06_scanpy_tradeoffs.py`. They use one immutable deterministic
CSR AnnData, the same observations, genes, grouping, expression source, 6 × 4.5
inch canvas, 80 DPI, and PNG output. The smoke test compares the prepared
embedding values, group means, and detection fractions before rendering all
eight matched figures.

## Equivalent examples

### Categorical embedding

```python
# ggann
ag.plot_embedding(
    adata, "umap", color="cell_type", pointdensity=False,
)

# Scanpy
sc.pl.embedding(
    adata, "umap", color="cell_type", show=False, return_fig=True,
)
```

### Gene-coloured embedding from a named layer

```python
# ggann: the selector makes the source local to this aesthetic
ag.plot_embedding(
    adata,
    "umap",
    color=ag.gene("NKG7", layer="logcounts"),
    pointdensity=False,
)

# Scanpy: the layer applies to the plotting call
sc.pl.embedding(
    adata,
    "umap",
    color="NKG7",
    layer="logcounts",
    use_raw=False,
    show=False,
    return_fig=True,
)
```

### Dotplot

```python
# ggann
ag.plot_dotplot(
    adata, genes, group_by="cell_type", use_raw=False,
    categories_order=cell_types,
)

# Scanpy
sc.pl.dotplot(
    adata, genes, "cell_type", use_raw=False,
    categories_order=cell_types, show=False, return_fig=True,
)
```

### Matrixplot

```python
# ggann
ag.plot_matrixplot(
    adata, genes, group_by="cell_type", use_raw=False,
    categories_order=cell_types,
)

# Scanpy
sc.pl.matrixplot(
    adata, genes, "cell_type", use_raw=False,
    categories_order=cell_types, show=False, return_fig=True,
)
```

Parity here means the same population and statistics, not pixel identity. The
two libraries deliberately use different plotting systems and defaults.

## Where Scanpy is stronger

Choose Scanpy when its canonical plot is sufficient and any of these dominate:

- minimum latency;
- constrained peak memory;
- an optimized direct Matplotlib implementation;
- familiar Scanpy defaults with little post-construction customization;
- large-sparse plotting speed rather than grammar composition.

The measured primary benchmark below favours Scanpy for preparation,
end-to-end time, and peak RSS. ggann therefore makes no "blazing fast" claim.

## Where ggann + annplyr is stronger

Choose ggann when the plot is a specification you expect to extend. For
example, this ordinary `plotnine.ggplot` adds a layer, facet, colour scale, and
theme without requiring a new ggann helper:

```python
(
    ag.gganndata(
        adata,
        aes(
            x=ag.obsm("umap", 0),
            y=ag.obsm("umap", 1),
            color=ag.gene("NKG7", layer="logcounts"),
            # Include the facet field in the projected plot data.
            group=ag.obs("condition"),
        ),
        # Two coordinates plus one gene; obs metadata is not charged.
        max_matrix_values=3 * adata.n_obs,
        add_theme=False,
    )
    + geom_point(size=1.8, alpha=0.85)
    + facet_wrap("condition")
    + scale_color_cmap(cmap_name="magma")
    + theme_classic()
)
```

This path is useful when you need:

- arbitrary plotnine layers, scales, annotations, facets, and themes;
- reusable publication figure specifications and multi-panel composition;
- explicit `obs()`, `gene()`, and `obsm()` resolution instead of precedence
  inferred from a bare string;
- different genes mapped from `.X`, `.raw`, or named layers in one grammar;
- only requested sparse or backed fields projected before conversion;
- a cumulative `max_matrix_values` rejection before extraction;
- positional safety for views, duplicate names, and reordered observations;
- annplyr tabular verbs and ggann plots over the same immutable object;
- a custom visual extension without waiting for a package-specific helper.

The materialization budget covers logical expression and `obsm` values, not
small `obs` metadata columns. annplyr typed errors propagate. Neither the
grammar call nor the plotting helpers mutate the input AnnData.

## Honest trade-off table

| Consideration | Standard Scanpy plotting | ggann + annplyr |
|---|---|---|
| Standard-plot speed | Currently faster in the matched primary suite | Slower; plotnine draw/save dominates |
| Peak memory | Lower in the measured large-CSR suite | Higher while plotnine and Matplotlib build artists |
| Composability | Function parameters and returned Matplotlib objects | Ordinary plotnine layers, scales, facets, guides, and themes |
| Source-selection clarity | Call-wide `layer=`/`use_raw=` plus gene keys | Bare-name convenience plus explicit `obs()`, `gene()`, and `obsm()` selectors |
| Sparse/backed extraction | Efficient canonical Scanpy paths | annplyr projects requested fields; sparse aggregation stays sparse |
| Materialization safeguards | No equivalent cumulative plotting argument | `gganndata(max_matrix_values=...)` and direct annplyr budgets |
| Custom themes and annotations | Matplotlib customization after or through plot parameters | Grammar-native and reusable in the plot specification |
| Return types | Vary by helper and `return_fig=` option | Ordinary `plotnine.ggplot` for grammar-native helpers; documented grid exceptions |
| Learning curve | Lowest for existing Scanpy users | Requires plotnine grammar concepts |
| Best fit | Fast canonical exploratory plots | Explicit extraction and deeply customized publication plots |

Neither path is universally superior. A practical workflow can use Scanpy for
fast standard diagnostics and ggann for the smaller set of figures whose
composition or source semantics justify the extra rendering cost.

## Measured benchmark

This table is generated from the versioned raw result
`benchmarks/results/scanpy-extended-csr.json`; it is not executed during the
documentation build.

```{include} ../_includes/scanpy-extended-csr.md
```

Regenerate the full seven-repeat comparison on a quiet machine with one
computational thread:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
MPLBACKEND=Agg MPLCONFIGDIR=/tmp/ggann-benchmark-mpl \
NUMBA_CACHE_DIR=/tmp/ggann-benchmark-numba PYTHONHASHSEED=0 \
PYTHONPATH=src python benchmarks/compare_scanpy.py \
  --preset extended --formats csr --workloads primary --sources x \
  --repeats 7 --seed 20260809 --include-cold-start \
  --isolated-memory-stages preparation,end_to_end \
  --isolated-memory-repeats 7 \
  --output benchmarks/results/scanpy-extended-csr.json \
  --report benchmarks/results/scanpy-extended-csr.md

python benchmarks/render_scanpy_vignette.py \
  benchmarks/results/scanpy-extended-csr.json \
  docs/_includes/scanpy-extended-csr.md
```

The small deterministic correctness/render smoke remains part of documentation
CI:

```bash
GGANN_DOCS_OFFLINE=1 MPLBACKEND=Agg PYTHONPATH=src \
  python examples/vignettes/06_scanpy_tradeoffs.py
```

## Complete executable source

```{literalinclude} ../../examples/vignettes/06_scanpy_tradeoffs.py
:language: python
:linenos:
```
