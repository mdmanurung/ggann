# Publication-quality multi-panel figure

Helper plots stay composable through final figure assembly. Build each panel,
use `compose` to arrange and tag them, then save once at the target physical
dimensions.

The executable example combines an embedding, dotplot, violin, and composition
summary. It writes to a temporary directory so the documentation build never
modifies committed figures.

## Executed source

```{literalinclude} ../../examples/vignettes/05_publication_panels.py
:language: python
:linenos:
```

For journal output, prefer PDF or SVG when the target accepts vectors. Set width
and height in `mm` to match the final layout, choose a single legible base font
size, and inspect the rendered output at its publication dimensions. Rasterize
only dense point layers when file size becomes impractical.
