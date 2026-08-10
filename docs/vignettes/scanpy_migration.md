# From Scanpy analysis to reusable figures

Suppose a Scanpy workflow has already produced a UMAP and curated PBMC labels.
The immediate review question is whether the annotation agrees with canonical
lineage markers. ggann keeps that familiar one-call plotting vocabulary, but
returns objects that can be styled, combined, inspected, and exported later.

The executable example uses the 2,638-cell `pbmc3k_processed` dataset -- the
PBMC sample behind Scanpy's and Seurat's clustering tutorials -- and six
canonical markers: CD3D, IL7R, MS4A1, NKG7, GNLY, and LYZ. CD3D, IL7R, and LYZ
did not survive highly-variable-gene selection, so the marker calls pass
`use_raw=True` to read all 13,714 genes from `.raw`.

## Familiar calls, reusable results

| Scanpy | ggann |
|---|---|
| `sc.pl.embedding(adata, "umap", color=...)` | `ag.plot_embedding(adata, "umap", color=...)` |
| `sc.pl.dotplot(adata, genes, groupby=...)` | `ag.plot_dotplot(adata, genes, group_by=...)` |
| `sc.pl.matrixplot(adata, genes, groupby=...)` | `ag.plot_matrixplot(adata, genes, group_by=...)` |
| `sc.pl.violin(adata, genes, groupby=...)` | `ag.plot_violin(adata, genes, group_by=...)` |

`group_by` is spelled consistently across ggann helpers. `show=` and
`return_fig=` are unnecessary because construction and rendering are separate:
every call below returns the plot without drawing it.

```python
sc.pl.embedding(
    adata,
    basis="umap",
    color="louvain",
    show=False,
)

embedding = ag.plot_embedding(
    adata,
    basis="umap",
    color="louvain",
    label=True,
)
```

Expression source choices are explicit in a real analysis. Here `.raw` contains
log-normalized expression, so both marker summaries request `use_raw=True`:

```python
genes = ["CD3D", "IL7R", "MS4A1", "NKG7", "GNLY", "LYZ"]

markers = ag.plot_dotplot(
    adata,
    genes,
    group_by="louvain",
    use_raw=True,
)
distribution = ag.plot_violin(
    adata,
    ["CD3D"],
    group_by="louvain",
    use_raw=True,
    add_box=True,
    stats=False,
)
```

The dotplot colour is the arithmetic mean of log-normalized expression and its
size is the fraction above zero. The violin retains every cell and adds a box
whose centre and interval are the median and interquartile range. No test is
requested in this descriptive review.

## Assemble the figure when the question is clear

Composition does not require extracting artists or redrawing individual
panels:

```python
figure = ag.compose(
    [embedding, markers],
    ncol=2,
    widths=(0.9, 1.3),
    gap=2,
    tag_levels=None,
)
```

`embedding`, `markers`, and `distribution` remain ordinary
`plotnine.ggplot` objects. Add a plotnine layer, scale, facet, label, coordinate
system, or theme to any one of them before composition. The executable workflow
also checks that the embedding contains all 2,638 cells, the marker table
contains all six genes, the violin contains all 2,638 measurements, and the
input AnnData
fingerprint is unchanged.

Visual equivalence with Scanpy means the same cells, variables, groups, and
statistics, not pixel identity between Matplotlib and plotnine. See
{doc}`../scanpy_parity` for the complete semantic mapping.

## Executed source

```{literalinclude} ../../examples/vignettes/01_scanpy_migration.py
:language: python
:linenos:
```
