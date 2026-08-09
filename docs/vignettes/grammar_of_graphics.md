# Ask a custom question with the grammar of graphics

The deterministic PBMC-like fixture contains control and stimulated cells. Its
stimulated cells have a higher MKI67 program by construction. The question is
whether that condition-associated expression pattern remains visible across
the embedding rather than being driven by one annotated lineage.

`gganndata` resolves only the fields referenced by `aes(...)`, then returns an
ordinary `plotnine.ggplot`. Explicit selectors make the source of every visual
variable reviewable:

```python
plot = gganndata(
    adata,
    aes(
        x=obsm("umap", 0),
        y=obsm("umap", 1),
        color=gene("MKI67", layer="logcounts"),
        group=obs("condition"),
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
3. `facet_wrap` separates control and stimulated cells;
4. `theme_classic` replaces the default theme.

`condition` is included as a `group` aesthetic so it is present in the resolved
table when the facet evaluates it. `add_theme=False` leaves the theme decision
to the final grammar. The executable check computes condition means from
`plot.data` and verifies that stimulated MKI67 is higher before rendering.

## Executed source

```{literalinclude} ../../examples/vignettes/02_grammar_of_graphics.py
:language: python
:linenos:
```

Use `plot.data` for transparent debugging or a methods audit: it is the exact
pandas table passed to plotnine. Treat it as returned plot state, not as a live
view of `AnnData`.
