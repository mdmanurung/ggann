# Statistical layers and pseudobulk

## Statistical layers

ggann re-exports statistical layers from
[`plotnine-extra`](https://github.com/mdmanurung/plotnine-extra). They compose
with plotnine-native ggann helpers.

```python
import scanpy as sc
import ggann as ag

adata = sc.datasets.pbmc3k_processed()
group = "louvain"

scatter = ag.plot_qc_scatter(
    adata,
    x="n_counts",
    y="n_genes",
) + ag.stat_cor()

violin = ag.plot_violin(
    adata,
    ["CD3D"],
    group,
    use_raw=True,   # CD3D lives in .raw, not in the highly-variable .X
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

## What each statistic shows

ggann draws these numbers; it does not choose or correct them. Every test
below is computed on the cells (or points) in the layer, so **cells are treated
as independent observations**. Cells from one donor are not independent, and
cell-level p-values are usually far too small for a claim about donors. For a
donor-level claim, aggregate with `pseudobulk` and test the profiles.

No multiple-testing correction is applied unless the table says so.

| Statistic | Test | Sides | p-value shown | Correction |
|---|---|---|---|---|
| `plot_violin(stats=True)`, `plot_box(stats=True)` | Two groups: Mann-Whitney U, unpaired. More than two: Kruskal-Wallis. One global test per feature facet. | Two-sided (Kruskal-Wallis has no direction) | Raw p, three significant digits. Computed on the cells that remain after `downsample=`. | None across facets |
| `stat_compare_means` | `method="wilcox.test"` (default), `"t.test"` (Welch), `"anova"`, `"kruskal.test"`. More than two groups with a pairwise method switch to Kruskal-Wallis or ANOVA. | Two-sided | Raw p, `p_digits=3` | None, including with `comparisons=` or `ref_group=` |
| `stat_pwc` | `method="wilcox.test"` (default) or `"t.test"`, all pairs | Two-sided | The default label `p.format` is the raw p. Use `label="p.adj.format"` or `"p.adj.signif"` for the adjusted p. | Yes: `p_adjust_method="holm"` by default; computed whether or not you show it |
| `stat_anova_test`, `stat_kruskal_test` | One-way ANOVA, Kruskal-Wallis | Not directional | Raw p | None |
| `stat_pvalue_manual` | None; it draws p-values you supply | As you computed | Your values | None |
| `stat_cor` | Pearson correlation by default (`method="spearman"` or `"kendall"`) | Two-sided (`alternative=`) | Raw p from the test, `p_accuracy=0.001`; needs at least 3 points | None |
| `stat_regline_equation` | Least-squares polynomial fit; no test | Not applicable | None. The label shows the equation and R-squared only. | Not applicable |
| `plot_volcano`, `plot_ma` | None; they threshold p-values you computed | Not applicable | `plot_volcano` plots `-log10` of `pvals_adj` from Scanpy's result. `plot_ma` thresholds the column named by `pval`. | Whatever produced the column |

The Wilcoxon rank-sum label is the Mann-Whitney U test from SciPy. SciPy picks
the exact p-value only when one sample has at most 8 values and there are no
ties. In every other case, which includes nearly all single-cell data, it uses
the normal approximation with a tie correction. R's `wilcox.test` uses the
exact p-value for samples below 50 without ties, so for samples of about 9 to 49
untied values the third digit of the p-value can differ from R's. With
hundreds of cells the two agree.

**Volcano and MA thresholds.** `plot_volcano(lfc=1.0, padj=0.05)` marks a gene
significant when `pvals_adj < padj` and `|logfoldchanges| >= lfc`.
`plot_ma(lfc=1.0, padj=0.05)` marks it significant when the `pval` column is
`< padj` and `|lfc_col| >= lfc`. Pass `lfc=0` to `plot_ma` to threshold on
`padj` alone. A gene with a missing p-value is never significant.

**Scanpy's marker tests.** `plot_volcano` and the rank-genes plots read the
result of `sc.tl.rank_genes_groups`. With `method="wilcoxon"` Scanpy tests each
group against the rest, two-sided, using the normal approximation without tie
correction by default. `pvals_adj` is Benjamini-Hochberg across genes within
each group (`corr_method`). The volcano therefore shows a gene-wise correction
and no correction across groups. Scanpy tests cells as independent, so these
p-values rank markers; they do not estimate replicate-level evidence.
`values_to_plot="pvals"` or `"pvals_adj"` in the rank-genes plots draws these
stored values unchanged.

For paired tests, prespecified contrasts, or a different label, use an explicit
`stat_compare_means(...)`, `stat_pwc(...)`, or `stat_pvalue_manual(...)` layer.

## Pseudobulk

`pseudobulk` aggregates cells to one profile per sample, optionally within each
group. It returns a new `AnnData` and requires the `pseudobulk` extra. Samples
or donors become the replicates: one row per `sample_col` value (and per
`group_by` level), so a test on the result compares donors, not cells. Cell-level
tests, by contrast, treat every cell as an independent replicate. `pseudobulk`
itself runs no test.

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

Pass `use_raw=True` instead of `layer=` to aggregate `counts_adata.raw`. The
counts and the gene set then both come from `.raw`, which may hold more genes
than `.var`, and the result carries the raw genes. `use_raw=True` without
`adata.raw` raises `ValueError`, as does combining it with `layer=`.
decoupler checks for integer counts unless `skip_checks=True` is set.

Profiles built from fewer than `min_cells` cells (default 10) are dropped, and a
`UserWarning` gives how many. Pass `min_cells=0` to keep them. The cell count of
each profile stays in `pb.obs["psbulk_cells"]`. Cells with a missing
`sample_col` or `group_by` value are dropped first, with their own warning.

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
