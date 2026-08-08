# Statistical layers and pseudobulk

## Statistical layers

ggann re-exports statistical layers from
[`plotnine-extra`](https://github.com/mdmanurung/plotnine-extra). They compose
with plotnine-native ggann helpers.

```python
import scanpy as sc
import ggann as ag

adata = sc.datasets.pbmc68k_reduced()
group = "bulk_labels"

scatter = ag.plot_qc_scatter(
    adata,
    x="n_counts",
    y="n_genes",
) + ag.stat_cor()

violin = ag.plot_violin(
    adata,
    ["CD3D"],
    group,
) + ag.stat_central_tendency()
```

| Layer | Purpose |
|---|---|
| `stat_compare_means`, `stat_pwc` | Group comparisons |
| `stat_pvalue_manual` | User-supplied p-values |
| `stat_cor`, `stat_regline_equation` | Correlation and regression labels |
| `stat_anova_test`, `stat_kruskal_test` | Omnibus tests |
| `stat_central_tendency` | Mean or median marker |
| `geom_signif` | Significance brackets |

`plot_violin` and `plot_box` accept `stats=True` as a shortcut for adding a
group-comparison layer.

## Pseudobulk

`pseudobulk` aggregates cells to one profile per sample, optionally within each
group. It returns a new `AnnData` and requires the `pseudobulk` extra.

The input below assumes that `counts_adata.obs` contains `donor` and
`cell_type`, and that `counts_adata.layers["counts"]` contains integer counts:

```python
pb = ag.pseudobulk(
    counts_adata,
    sample_col="donor",
    group_by="cell_type",
    layer="counts",
    mode="sum",
)
```

Pass `use_raw=True` instead of `layer=` to aggregate `counts_adata.raw`.
`raw=` is accepted as a deprecated spelling. decoupler checks for integer
counts unless `skip_checks=True` is set. Profiles with fewer than `min_cells`
cells are omitted.

The returned object works with the grammar and expression-summary helpers when
their required observation columns are present:

```python
from ggann import aes, gene, gganndata
from plotnine import geom_boxplot

boxplot = (
    gganndata(pb, aes("cell_type", gene("CD3D"), fill="donor"))
    + geom_boxplot()
)

dotplot = ag.plot_dotplot(
    pb,
    ["CD3D", "NKG7", "CST3"],
    group_by="cell_type",
    use_raw=False,
)
```
