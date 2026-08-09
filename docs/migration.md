# Migration guide

## Release-candidate change ledger

No breaking change to an already released ggann public API is recorded for the
current `0.1.0` candidate. The release-candidate work adds optional grammar and
coordinate-extraction budgets. Before publishing another version, update this
table from the actual API diff and add a deprecation path for every renamed or
removed symbol.

| Before | Current | Status | Action |
|---|---|---|---|
| `gganndata(...)` without a hard read budget | `gganndata(..., max_matrix_values=None)` | Additive | Set an integer only when the mapping needs a hard cumulative limit |
| `embedding_coords(adata, basis, n=2)` | Adds keyword-only `max_matrix_values=None` | Additive | Existing calls are unchanged |
| Grouped summaries could retain a missing (`NA`) grouping row | Rows missing any `group_by`/`split_by` key are omitted | Consistency fix | Fill missing labels before plotting if they should form an explicit category |
| Sparse grouped means had backend-dependent `NaN` behavior | Expression `NaN` values are skipped for means and treated as not detected for fractions | Correctness fix | No change for finite matrices |

Documentation-only clarification is not a public API change. The current
contracts are summarized below so later releases have an explicit baseline.

## Current naming baseline

| Concept | Canonical ggann spelling | Notes |
|---|---|---|
| Primary group | `group_by` | Scanpy spells this `groupby` |
| Secondary split | `split_by` | Facet or within-group split, depending on helper |
| Aesthetic colour | `color` | `scale_colour_*` aliases remain available |
| Expression layer | `layer` | Mutually exclusive with `use_raw=True` |
| Raw matrix switch | `use_raw` | `None` selects `.raw` when present |
| Cell cap | `downsample` | Never applied silently |
| Sampling seed | `random_state` | Defaults to deterministic `0` |
| Primary category order | `categories_order` | Must include observed non-missing groups |
| Grammar read budget | `max_matrix_values` | Cumulative expression and `obsm` values; observation metadata is free |

## Moving from Scanpy

The closest direct translations are:

```python
# Scanpy
sc.pl.dotplot(adata, genes, groupby="cell_type", show=False)

# ggann
plot = ag.plot_dotplot(adata, genes, group_by="cell_type")
```

ggann returns rather than immediately showing a plot. Compose it with plotnine,
then draw or save explicitly:

```python
from plotnine import labs, theme

plot = plot + labs(title="Markers") + theme(figure_size=(7, 4))
plot.save("markers.pdf", width=180, height=100, units="mm")
```

See {doc}`vignettes/scanpy_migration` and {doc}`scanpy_parity` for full mappings
and semantic differences.

## Source-selection baseline

With neither `layer` nor `use_raw` set, ggann uses `.raw` when available and
otherwise `.X`. Use `use_raw=False` to force `.X`. A named layer always wins;
combining `layer=` with `use_raw=True` is an error.

Grammar code that previously relied on a bare-name collision should move to an
explicit selector:

```python
from ggann import gene, obs

cell_annotation = obs("CD3D")
gene_expression = gene("CD3D", use_raw=True)
```

## Return-type baseline

All plotnine-native helpers return `plotnine.ggplot`. The only plotting
exceptions are `plot_clustermap` and `plot_upset`, which return their grid
backend objects. Code should not assume those two support plotnine's `+`
operator.

## Deprecation policy for future changes

An unavoidable public rename should retain the old spelling for at least one
minor release, emit a targeted `DeprecationWarning`, document the replacement,
and add both old and new calls to the API contract tests. Silent behavior changes
are not an acceptable migration strategy.
