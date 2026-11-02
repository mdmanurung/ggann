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
sans-serif family.

If Arial and Helvetica are missing, ggann warns once per family and process,
and Matplotlib uses DejaVu Sans. SVG text stays editable and still names the
requested family. PDFs embed the font that was actually used. Text widths
differ between fonts, so check the final layout on the machine that exports.
An unknown family in `set_theme(family=...)` warns the same way and does not
raise. Only a simulated Arial was tested (DejaVu Sans registered under the
name Arial), because the test host has no Arial.

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
normal exit or an exception. A `set_theme` call or an `rcParams` edit made inside
the block is discarded on exit. Plots keep the style specification captured when
they were constructed, so export can happen outside the block.

The active style is not inherited by new threads. A new thread sees the
process-wide plotnine theme and `rcParams` but not the style. Concurrent use of
`style_context` from several threads is unverified, so run publication exports
from one thread.

`set_theme` affects plots constructed afterwards. `reset_theme` installs
`theme_gray` and resets `ggann.sizes` to a base size of 11. It does not restore
an earlier custom theme; use `style_context` for that. `set_theme`,
`theme_ggann` and `sizes.update` raise `ValueError` for a `base_size` that is
not a finite positive number, before they change any state.

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
from zero.

Only those first eight colours are measured. A larger palette extends a fixed
colour sequence, so a 21st category never recolours the first 20. Separation
beyond eight is not guaranteed. Measured deuteranopia minimum separation is 8.6
at nine colours, 4.0 at 10 to 14 and 3.5 at 15 to 20. ggann warns once per
palette when stored or generated colours fall below 5 CIEDE2000 under any
simulation, which starts at 10 categories on the default palette. It also warns
once when stored colours are duplicated. The warning only reports: colours are
kept as given. Add direct labels, shapes, or faceting for those figures.

Colours follow category order. Reorder the categories and the colours move with
them.

Outside publication mode, a discrete scale without stored
`adata.uns["<col>_colors"]` uses Scanpy's default palette (20, 28 or 102
colours), so colours match `sc.pl.umap`. More than 102 categories use plotnine
hues. Stored colours that are invalid warn and fall back to the default palette.
In the raw grammar, `gganndata` adds no default colour scale: add
`+ scale_color_obs(adata, "cell_type")` to use the stored palette. Helpers apply
it automatically.

## Colour maps for expression

Expression helpers take `cmap=None`. Precedence is an explicit `cmap`, then the
style's sequential map inside `style_context`, then the helper's own default
(`magma` for `plot_features` and `plot_density`, `Reds` for dotplots, `viridis`
for the rest). Pass a name to override the style, including `"Reds"`. With
`standard_scale="zscore"`, an explicit `cmap` is drawn with zero-centred limits.
In the Matplotlib backend an explicit `cmap` is drawn without the zero-centred
norm. In `plot_embedding`, `low` and `high` only drive the legacy two-colour
gradient and are ignored inside a style context.

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
`"i"` sequences are available. Use `compose(tag_levels=)` for panel letters:
plotnine-extra's `plot_annotation` has no tag option. `compose` and `tag_panels`
raise `TypeError` for anything that is not a plot.

`guides="collect"` removes only later guides whose mapping, scale class, limits,
breaks, labels, missing-value handling, and palette match exactly. It does not
merge merely similar biological encodings. The remaining legend stays inside its
panel; it is not moved to the figure margin.

Save a composition with `.save(path, width=, height=, units=, dpi=)`. Any
positive integer `dpi` works, including 150; the default is 300. Plotnine's
`ggsave` does not accept compositions.

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

`height` is required. The `dpi` argument takes any positive integer for PNG and
TIFF (150 works for drafts); bool, zero, negative, non-finite and non-integer
values raise. Styles still accept only 300 or 600 in `PublicationStyle(dpi=)`.
Raster size is `floor(size_in_inches × dpi)` pixels, so 89 mm at 150 dpi is 525
px. SVG is written first when requested and uses text nodes; PDF requests
embedded Type 42 fonts. Backgrounds may be white, transparent, or any valid
Matplotlib colour. Temporary figures are closed, and repeated exports do not
mutate the plot object.

Only svg, pdf, png, tiff and tif count as a suffix. Any other dot stays in the
stem, so `fig1.v2` writes `fig1.v2.svg`, `fig1.v2.pdf` and `fig1.v2.png`. A
suffix that is missing from an explicit `formats` raises: `x.png` with
`formats=("svg",)` is an error. `.jpg` and `.eps` raise.

`save_publication` warns when text, a legend or a colourbar extends beyond the
final canvas or is clipped by its axes, and names the items and the overshoot in
millimetres. It never crops, resizes or enlarges the canvas, so fix the figure
and export again. The check adds about 0.1 s per call.

The same exporter accepts the objects returned by `plot_clustermap` and
`plot_upset`, even though those two documented exceptions are not plotnine
objects. Both export to SVG, PDF and PNG.

A clustermap is fitted to the canvas at export. PyComplexHeatmap lays out the
body, dendrograms, annotations and legends for the figure size at construction
and leaves no room for tick labels or its own legends, so `save_publication`
measures them and shrinks the heatmap body until every label, legend and
colourbar is inside. Dendrograms, annotations and legends keep their size and
follow the body. The canvas and dpi are not changed. This also applies to
`TidyHeatmap.save`.

Fitting stops when the body would be smaller than 20 mm. The overflow warning
then stays, and says the heatmap was shrunk to its minimum body size. On
`pbmc68k_reduced` with six genes, 183 × 120 mm now fits with nothing cut off
(it cut off 39 items, by up to 32 mm, before). At 89 × 70 mm the ten long group
names and the two-column legend do not fit, so it warns. Shorten the labels or
pass `plot_legend=False`.

Construct the clustermap at the export size when you can, so the proportions of
the annotations and dendrograms start right:

```python
import matplotlib.pyplot as plt

plt.figure(figsize=(89 / 25.4, 70 / 25.4))
heatmap = ag.plot_clustermap(adata, genes, group_by="cell_type")
ag.save_publication(heatmap, "clustermap", width="single-column", height=70)
```

A clustermap built at its own size and resized at export has distorted
annotation and dendrogram proportions, even though nothing is cut off.

A Matplotlib dotplot or matrixplot measures its left margin from the rendered
y labels, at export size, and its legend fits the 89 mm column on
`pbmc68k_reduced`. The plotnine backend can still report a legend overshoot of a
fraction of a millimetre at 89 mm with long labels.

## Rasterize data, not the annotation system

Dense scatter and tile helpers expose `rasterized=False`. Set it to true only
when the vector data layer makes a file impractically large. An SVG made with
`rasterized=True` holds one embedded `<image>` for the data and vector `<text>`
for every label, in both backends; without it the SVG has no image:

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

Only the point or tile artists become raster data. Text, axes, legends,
colourbars, tags, and cell annotations remain editable vectors. `annotate="auto"` labels a
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
