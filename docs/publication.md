# Publication figures

`ggann` separates scientific preparation from figure-wide styling. The same
helper call can serve exploration and final assembly; publication mode changes
themes, palettes, guides, axes, and rendering defaults without changing the
prepared table or dropping observations.

The generic presets are journal-oriented starting points, not promises of
compliance with a named venue. Confirm the destination journal's current width,
font, and file-format requirements before submission.

## Define the figure contract first

Before styling, write down four things:

1. The scientific claim the figure supports.
2. A panel map in which each panel has one job.
3. The population and `n` represented in each panel.
4. Every centre, interval, model, and statistical test shown.

The {doc}`vignettes/publication_panels` workflow applies that contract to the
bundled PBMC68k data. Its manifest records `n=700`, arithmetic means and
detection fractions for the dotplot, median and interquartile range for the
violin, descriptive-only panels with no inferential tests, the exact canvas,
the shared colour vocabulary, and a SHA-256 input fingerprint.

## Styles and presets

`PublicationStyle` is a frozen, validated dataclass. Build one from a generic
width preset and override only what the figure needs:

```python
import ggann as ag

single = ag.publication_style("single-column")  # 89 mm, default height 70 mm
double = ag.publication_style("double-column")  # 183 mm, default height 120 mm
custom = ag.publication_style(
    double,
    base_size=7,
    tag_size=8.5,
    line_width=0.55,
)
```

Text defaults span 5.5–7 points, tags are 8 points, and primary lines are 0.5
points. Font lookup follows Arial, Helvetica, DejaVu Sans, then the generic
sans-serif family. This makes the fallback explicit while allowing the figure
to render on systems without proprietary fonts.

Use `style.to_dict()` when a pipeline needs to write the complete style into a
manifest.

## One theme, four axis treatments

`theme_publication` is a normal plotnine theme. Add it to an existing plot, and
add any plot-specific override after it:

```python
from plotnine import theme

final = exploratory + ag.theme_publication(axes="embedding")
final = final + theme(legend_position="bottom")
```

| `axes=` | Intended use | Treatment |
|---|---|---|
| `"standard"` | bars, trends, QC, differential expression | restrained classic axes |
| `"embedding"` | UMAP, t-SNE, PCA | hides arbitrary numeric units and axis lines |
| `"matrix"` | dotplots, heatmaps, correlations | removes redundant axes and anchors rotated labels |
| `"distribution"` | violin, box, sina, ridge | rotates category labels for compact layouts |

A later user theme, scale, coordinate system, layer, facet, label, or annotation
remains authoritative.

## Apply one style across a workflow

`style_context` coordinates helper defaults, bare plotnine plots, the optional
direct Matplotlib backend, clustermap and UpSet specifications, global plotnine
theme state, `ggann.sizes`, and Matplotlib `rcParams`:

```python
with ag.style_context(custom):
    umap = ag.plot_embedding(adata, "umap", color="cell_type")
    markers = ag.plot_dotplot(adata, genes, group_by="cell_type")
    figure = ag.compose([umap, markers], ncol=2, widths=(0.9, 1.3), gap=2)
```

Contexts can be nested. The current plotnine theme, every `ggann.sizes` field,
all Matplotlib parameters, and the enclosing context are restored exactly after
normal exit or an exception. Plots keep the style specification captured when
they were constructed, so export can happen outside the block.

Outside a context, helpers retain their standard exploratory appearance. This
makes publication styling opt-in and keeps notebooks concise.

## Stable, validated colour vocabularies

Use `publication_palette` for qualitative, sequential, or zero-centred
diverging values:

```python
cell_types = ["T cell", "B cell", "NK cell", "Monocyte"]
cell_type_colours = ag.publication_palette(
    "qualitative", categories=cell_types
)
expression_colours = ag.publication_palette("sequential", n=64)
effect_colours = ag.publication_palette("diverging", n=65)
```

Supplying categories returns an insertion-ordered category-to-hex mapping. Use
that mapping once in `adata.uns["cell_type_colors"]` when the vocabulary should
travel with the annotation.

In publication mode, a stored AnnData palette is reused only if it has exactly
one syntactically valid colour for every categorical level. A missing or invalid
palette produces a deterministic fallback; invalid input also emits a targeted
warning. `ggann` never repairs or mutates `adata.uns` implicitly. Missing
observations remain represented with neutral grey `#B3B3B3`.

The core eight qualitative colours are checked in normal vision, grayscale,
protanopia, deuteranopia, and tritanopia simulations. The test requires minimum
CIEDE2000 separation of 10 normally and 5 in each colour-vision simulation.
Sequential luminance is monotonic, and both diverging arms are monotonic away
from zero. Vocabularies beyond eight categories need figure-specific review or
a redundant encoding such as shape, direct labels, or faceting.

## Weighted layout, gaps, guides, and tags

`compose` accepts positive relative column widths and row heights plus a
physical gap in millimetres:

```python
figure = ag.compose(
    [umap, dotplot, violin, proportions],
    ncol=2,
    widths=(0.95, 1.25),
    heights=(1.0, 1.0),
    gap=2,
    guides="keep",
    tag_levels="auto",
)
```

The return value remains an instance of plotnine's `Compose`. In publication
mode, automatic tags are bold lowercase; explicit `"A"`, `"a"`, `"1"`, and
`"i"` sequences are available. `guides="collect"` removes only later guides
whose mapping, scale class, limits, breaks, labels, missing-value handling, and
palette match exactly. It does not merge merely similar biological encodings.

## Exact, editable export

`save_publication` draws a copied object in a local font/rendering context and
never tight-crops the canvas:

```python
paths = ag.save_publication(
    figure,
    "figure_1",
    width="double-column",
    height=120,
    formats=("svg", "pdf", "png", "tiff"),
    dpi=600,
    background="white",
)
```

| Request | Result |
|---|---|
| suffixless stem, no `formats` | SVG, PDF, and PNG |
| suffixed filename, no `formats` | that one SVG, PDF, PNG, or TIFF |
| explicit `formats` | requested formats, emitted in SVG/PDF/PNG/TIFF order |
| numeric width | interpreted in `units` |
| `"single-column"` or `"double-column"` width | exact 89 or 183 mm preset |

Raster output accepts 300 or 600 DPI. SVG is written first when requested and
uses text nodes; PDF requests embedded Type 42 fonts; PNG and TIFF dimensions
are `round(size / 25.4 × dpi)` for millimetres. Backgrounds may be white,
transparent, or any valid Matplotlib colour. Temporary figures are closed, and
repeated exports do not mutate the plot object.

The same exporter accepts the objects returned by `plot_clustermap` and
`plot_upset`, even though those two documented exceptions are not plotnine
objects.

## Rasterize data, not the annotation system

Dense scatter and tile helpers expose `rasterized=False`. Set it to true only
when the vector data layer makes a file impractically large:

```python
with ag.style_context("double-column"):
    umap = ag.plot_embedding(
        adata,
        "umap",
        color="cell_type",
        rasterized=True,
    )
    matrix = ag.plot_matrixplot(
        adata,
        genes,
        group_by="cell_type",
        rasterized=True,
        annotate="auto",
    )
```

Only point or tile artists become raster data. Text, axes, legends, colourbars,
tags, and cell annotations remain editable vectors. `annotate="auto"` labels a
matrix cell only when it is at least 12 points wide and high at the final output
size; the label switches between black and white according to fill contrast.
Use `annotate="force"` when every value must appear, then inspect overlaps at
the final physical size.

## Final inspection checklist

- View the SVG or PDF at the submitted dimensions, not just zoomed on screen.
- Confirm every panel represents the intended population and report its `n`.
- Define all centres, intervals, detection thresholds, and tests in the legend.
- Check every text, tag, legend, and colourbar is inside the canvas.
- Review grayscale and colour-vision simulations; do not rely on colour alone.
- Verify vector text remains text and any rasterization is limited to data layers.
- Keep a machine-readable manifest with dimensions, DPI, style, palette, input
  fingerprint, software versions, and output hashes.

The centralized sizing, restrained typography, and editable-export principles
are informed by [CNSPlots](https://github.com/faridrashidi/cnsplots), without a
runtime dependency or copied implementation. The claim-first panel map and
final-size inspection follow the practical workflow used for high-impact
journal figure preparation.
