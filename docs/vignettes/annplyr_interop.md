# Plot annplyr pipelines with ggann

An annplyr result goes straight into ggann. `ap.tidy(adata)` returns a
`TidyAnnData`, an `AnnData` subclass whose verbs are methods, so every verb that
returns an object returns something `gganndata` and the `plot_*` helpers accept.
Verbs that return a table (`to_wide`, `join_features`, `summarize`, `count`)
return a DataFrame that plotnine draws directly. No adapter sits between them.

Names carry their container: `obs.louvain`, `raw.MS4A1`, `X.CST3`, `X_umap`.
Import them from `annplyr.data`.

```python
import annplyr as ap
from annplyr.data import X, X_umap, obs, raw

import ggann as ag

adata = ap.tidy(adata)
```

## Filter, derive, plot

A column made with `mutate` is an `obs` column, so it works as `color=` in a
helper and in `aes` in the grammar. Expression is read from `.raw` once, by
annplyr; ggann never sees the gene columns.

```python
scored = adata.filter(obs.percent_mito < 0.03).mutate(
    b_score=(raw.MS4A1 + raw.CD79A) / 2,
)
scored.pipe(ag.plot_embedding, color="b_score")
ag.gganndata(scored, aes("UMAP_1", "UMAP_2", color="b_score")) + geom_point()
```

`plot_embedding` sorts cells by value so high scores draw on top. Align its data
to the object by cell name, not by position. `rename`, `slice_max(..., by=)` and
`aggregate_cells` results feed ggann the same way.

## Embedding and gene tables for plain plotnine

`to_wide` names columns the way ggann resolves them: genes keep their bare
names and `X_umap` becomes `UMAP_1`, `UMAP_2`. Columns come out in source order
(`obs`, then matrices, then `obsm`), not in the order you list them.

```python
table = adata.to_wide(obs.louvain, X_umap, X.CST3, X.NKG7, max_matrix_values=adata.n_obs * 4)
ggplot(table, aes("UMAP_1", "UMAP_2", color="louvain")) + geom_point() + ag.theme_ggann()
```

`join_features(..., shape="long")` gives one row per cell and feature, with
`.cell`, `.feature` and `.abundance_raw`. That column has a pandas sparse dtype
and plotnine cannot compute box statistics on it (`ValueError: Cannot get a null
value for type: <class 'numpy.float32'>`). Densify the one column first:

```python
long = adata.join_features(raw.CST3, raw.NKG7, shape="long")
long[".abundance_raw"] = long[".abundance_raw"].sparse.to_dense()
ggplot(long, aes("louvain", ".abundance_raw", fill=".feature")) + geom_boxplot()
```

`join_features(..., shape="wide")` instead puts the genes in `obs`. If a gene
shares its name with an `obs` column, ggann uses the column and warns; write
`gene("CST3")` to force the expression.

## A custom depth-aware summary

Library-size normalization does not fully remove the sequencing depth effect.
`depth` is a median split of `n_counts`. The summary compares mean CST3 between
shallow and deep cells within each PBMC3K cluster. In CD14+ monocytes, the
largest myeloid cluster, mean CST3 is 2.88 in deeply sequenced cells against
2.06 in shallow ones. The script asserts on that cluster only, because the
rarest clusters hold as few as one cell per depth bin.

```python
summary = adata.summarize(
    mean_CST3=ap.mean(X.CST3),
    by=(obs.louvain, obs.depth),
    max_matrix_values=adata.n_obs,
)
```

A grouped expression summary can carry a pandas sparse dtype too. The script
densifies it before `geom_col` when it does.

## Counts, palettes and category order

`count(obs.a, obs.b)` lists groups in key order with an `n` column, and
`plot_proportions(group_by=..., split_by=...)` reports the same numbers; the
script merges the two and asserts they match.

`ap.store_palette(adata, "louvain", {...})` writes `uns["louvain_colors"]`.
`scale_color_obs`, `obs_colors` and `plot_embedding` read it, so the drawn
colours are the stored ones. `filter` drops unused categories, keeps the
remaining levels in their original order, and trims the stored palette to
match. A label keeps its colour after filtering.

## Executed source

The script below runs every snippet on PBMC3K, checks the plotted numbers
against the source object, and verifies that neither annplyr nor ggann edits the
input.

```{literalinclude} ../../examples/vignettes/04_annplyr_interop.py
:language: python
:linenos:
```

Select features explicitly on large objects. Pass `max_matrix_values=` when
exceeding a known memory boundary must be an error rather than a performance
surprise; `gganndata(..., max_matrix_values=...)` enforces the same budget.
