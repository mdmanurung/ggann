# Ask a custom question with the grammar of graphics

`NKG7` marks the cytotoxic granule program of NK and CD8 T cells. The question
is whether that program stays confined to the lymphoid side of the PBMC3K
embedding, or also appears among myeloid cells. `compartment` splits the
published Louvain labels into their lymphoid and myeloid lineages, so the facet
asks the biological question directly.

`gganndata` resolves only the fields referenced by `aes(...)`, then returns an
ordinary `plotnine.ggplot`. Explicit selectors make the source of every visual
variable reviewable:

```python
plot = gganndata(
    adata,
    aes(
        x=obsm("umap", 0),
        y=obsm("umap", 1),
        color=gene("NKG7", layer="logcounts"),
        group=obs("compartment"),
    ),
    max_matrix_values=3 * adata.n_obs,
    add_theme=False,
)
```

The materialization budget is exactly the two embedding coordinates plus one
gene for every observation. Observation metadata does not count toward the
matrix-value budget.

The workflow below demonstrates the four composition points expected from a
grammar interface:

1. `geom_point` adds a layer;
2. `scale_color_cmap` controls the continuous expression scale;
3. `facet_wrap` separates lymphoid and myeloid cells;
4. `theme_classic` replaces the default theme.

`compartment` is mapped to the non-visual `group` aesthetic, which puts it in the
resolved table from the start. A facet that names a field the first mapping
omitted would also work: `gganndata` reads late fields when they are added, and
counts them against `max_matrix_values`. `add_theme=False` leaves the theme
decision to the final grammar. The plot has no stored-colour scale, because
`gganndata` adds none; add `+ scale_color_obs(adata, key)` for a categorical
column with stored colours. The executable check computes compartment means
from `plot.data` and confirms the lymphoid excess (0.63 versus 0.20 mean
log-normalized `NKG7`) before rendering. The figure shows where that excess
sits: one bright island of NK and CD8 T cells, not a diffuse lymphoid shift.

## Executed source

```{literalinclude} ../../examples/vignettes/02_grammar_of_graphics.py
:language: python
:linenos:
```

Use `plot.data` for transparent debugging or a methods audit: it is the exact
pandas table passed to plotnine. Treat it as returned plot state, not as a live
view of `AnnData`.
