# Grammar-of-graphics workflow

`gganndata` resolves only the fields referenced by `aes(...)`, then returns a
plain `plotnine.ggplot`. Everything after that point is ordinary plotnine.

The workflow below demonstrates the four composition points expected from a
grammar interface:

1. `geom_point` adds a layer;
2. `scale_color_brewer` controls a scale;
3. `facet_wrap` creates panels;
4. `theme_classic` replaces the default theme.

`condition` is included as a `group` aesthetic so it is present in the resolved
table when the facet evaluates it. `add_theme=False` prevents ggann's default
theme from being added before `theme_classic`. The two embedding coordinates
are bounded by `max_matrix_values=2 * adata.n_obs`; the two observation fields
do not count toward that annplyr-compatible budget.

## Executed source

```{literalinclude} ../../examples/vignettes/02_grammar_of_graphics.py
:language: python
:linenos:
```

Use `plot.data` for transparent debugging: it is the exact pandas table passed
to plotnine. Treat it as returned plot state, not as a live view of `AnnData`.
