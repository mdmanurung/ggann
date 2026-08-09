# Quickstart

This guide uses Scanpy's bundled 700-cell PBMC68k subset. It has a UMAP,
curated `bulk_labels`, cell-cycle `phase`, scaled values in `.X`, and
log-normalized expression in `.raw`, so the core PBMC workflow runs as written
without downloading data.

```python
import scanpy as sc
import ggann as ag

adata = sc.datasets.pbmc68k_reduced()
genes = ["CD3D", "MS4A1", "NKG7", "GNLY", "CST3"]
```

## One-call biological summaries

Use a helper when the intended summary is standard:

```python
embedding = ag.plot_embedding(
    adata,
    basis="umap",
    color="bulk_labels",
    label=True,
)

dotplot = ag.plot_dotplot(
    adata,
    genes,
    group_by="bulk_labels",
    use_raw=True,
)

violin = ag.plot_violin(
    adata,
    ["CD3D"],
    group_by="bulk_labels",
    use_raw=True,
    add_box=True,
)
```

The embedding retains all 700 observations. The dotplot colour is arithmetic
mean log-normalized expression from `.raw`; dot size is the fraction above the
expression cutoff. The violin uses every `.raw` observation; its inner box
shows the median and interquartile range. These calls return ordinary
`plotnine.ggplot` objects.

## Compose unequal panels

`compose` adds physical gaps, relative row/column sizes, and panel tags while
remaining a native plotnine composition:

```python
figure = ag.compose(
    [embedding, dotplot],
    ncol=2,
    widths=(0.9, 1.3),
    gap=2,                 # millimetres
    tag_levels="auto",
    guides="keep",
)
```

`tag_levels="auto"` gives exploratory figures uppercase tags and publication
figures bold lowercase tags. `guides="collect"` suppresses a later guide only
when its mapping, limits, labels, and palette exactly match an earlier guide;
otherwise both remain visible.

## Switch the same calls to publication mode

Build panels inside `style_context`. The context applies final-size typography,
line weights, palette validation, family-specific axes, Matplotlib font
settings, and direct-backend defaults. It restores all global state on exit,
including after an exception.

```python
with ag.style_context("double-column"):
    embedding = ag.plot_embedding(
        adata, "umap", color="bulk_labels", label=True
    )
    dotplot = ag.plot_dotplot(
        adata, genes, group_by="bulk_labels", use_raw=True
    )
    figure = ag.compose(
        [embedding, dotplot],
        ncol=2,
        widths=(0.9, 1.3),
        gap=2,
    )

outputs = ag.save_publication(
    figure,
    "pbmc_markers",        # suffixless: SVG, PDF, and PNG
    width="double-column", # 183 mm
    height=90,
    dpi=600,
)
```

`outputs` is an ordered tuple of paths. A suffixed path such as
`"pbmc_markers.svg"` writes one format; `formats=("svg", "pdf", "png",
"tiff")` requests an explicit set. PNG and TIFF accept 300 or 600 DPI. The
canvas is never tight-cropped, so 183 × 90 mm at 600 DPI is reproducible.

Dense point or tile layers can be rasterized without flattening the rest of a
vector figure:

```python
with ag.style_context():
    dense = ag.plot_embedding(
        adata, "umap", color="bulk_labels", rasterized=True
    )

ag.save_publication(dense, "dense_umap.svg", height=70)
```

The SVG contains an image only for the point layer; text, axes, guides, and
annotations remain vectors. Leave `rasterized=False` for fully vector output.

## Use a stable colour vocabulary

Helpers reuse `adata.uns["<obs>_colors"]` only when it contains exactly one
valid colour per category. In publication mode an invalid stored palette emits
a warning and falls back deterministically without modifying `adata`.

For a custom annotation, make the vocabulary explicit once:

```python
import pandas as pd

cell_types = ["T cell", "B cell", "NK cell", "Monocyte", "Dendritic", "CD34+"]
colours = ag.publication_palette("qualitative", categories=cell_types)

adata = adata.copy()
adata.obs["cell_type"] = adata.obs["bulk_labels"].map(
    {
        "CD4+/CD25 T Reg": "T cell",
        "CD4+/CD45RA+/CD25- Naive T": "T cell",
        "CD4+/CD45RO+ Memory": "T cell",
        "CD8+ Cytotoxic T": "T cell",
        "CD8+/CD45RA+ Naive Cytotoxic": "T cell",
        "CD14+ Monocyte": "Monocyte",
        "CD19+ B": "B cell",
        "CD34+": "CD34+",
        "CD56+ NK": "NK cell",
        "Dendritic": "Dendritic",
    }
)
adata.obs["cell_type"] = adata.obs["cell_type"].astype(
    pd.CategoricalDtype(cell_types, ordered=True)
)
adata.uns["cell_type_colors"] = [colours[value] for value in cell_types]
```

Missing observations remain represented in embedding plots with the neutral
publication grey. The reviewed qualitative core contains eight colours; larger
vocabularies require visual review or a redundant encoding such as shape or
direct labels.

## Build a custom grammar plot

Use `gganndata` when plotnine already expresses the figure clearly. Include a
field in the aesthetic mapping when a later facet or layer needs that column;
`group=` is a convenient non-visual mapping for `phase` here.

```python
from ggann import aes, gganndata, obs, obsm
from plotnine import facet_wrap, geom_point, theme

plot = (
    gganndata(
        adata,
        aes(
            x=obsm("umap", 0),
            y=obsm("umap", 1),
            color=obs("bulk_labels"),
            group=obs("phase"),
        ),
    )
    + geom_point(size=1.0, alpha=0.8)
    + facet_wrap("phase")
    + theme(legend_position="bottom")
)
```

The returned object accepts any compatible plotnine geom, stat, scale,
coordinate system, facet, label, or theme. A bare plotnine plot built inside
`style_context` receives the publication current theme too.

## Choose the expression source

The same source arguments apply to the grammar and expression helpers:

```python
# scaled values in adata.X
x_plot = ag.plot_embedding(adata, "umap", color="CD3D", use_raw=False)

# log-normalized values in adata.raw.X
raw_plot = ag.plot_embedding(adata, "umap", color="CD3D", use_raw=True)

# a named layer on an object that provides one
counts_plot = ag.plot_embedding(your_adata, "umap", color="CD3D", layer="counts")
```

With neither argument, expression uses `adata.raw` when present and otherwise
`.X`, matching Scanpy's plotting convention. `layer=` and `use_raw=True` are
mutually exclusive. Pin individual genes with `ag.gene(...)` when one grammar
plot intentionally mixes sources.

## Bound grammar materialization

`gganndata` can reject a mapping before its first matrix read when the complete
request exceeds a known boundary. The following request contains two `obsm`
coordinates and one gene, or `3 × n_obs` logical matrix values; observation
metadata is not charged.

```python
bounded = gganndata(
    adata,
    aes(
        x=obsm("umap", 0),
        y=obsm("umap", 1),
        color=ag.gene("CD3D", use_raw=False),
    ),
    max_matrix_values=3 * adata.n_obs,
) + geom_point()
```

An invalid or exceeded budget raises `annplyr.AnnplyrError`. Continue with
{doc}`publication` for the full style, palette, rasterization, and editable
export contracts, or {doc}`concepts` for resolution, ordering, ownership, and
missing-value behavior.
