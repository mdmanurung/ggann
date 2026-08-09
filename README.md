# ggann

`ggann` turns `AnnData` into composable plotnine graphics and concise
single-cell plotting helpers. Start with a one-call biological summary, refine
it with the grammar of graphics, and export an editable multi-panel figure at
its final physical size.

[Documentation](https://mdmanurung.github.io/ggann/) ·
[Quickstart](https://mdmanurung.github.io/ggann/quickstart.html) ·
[API concepts](https://mdmanurung.github.io/ggann/concepts.html) ·
[API reference](https://mdmanurung.github.io/ggann/api.html)

## First plot

Given an `AnnData` object named `adata`, a useful plot is one call away:

```python
import ggann as ag

plot = ag.plot_embedding(adata, "umap", color="cell_type")
plot
```

Helpers prepare common single-cell summaries. The grammar entry point returns a
regular `plotnine.ggplot`, so plotnine layers, scales, facets, and themes work
without adapters:

```python
from ggann import aes, gganndata
from plotnine import facet_wrap, geom_point, theme_classic

plot = (
    gganndata(adata, aes("UMAP_1", "UMAP_2", color="cell_type", group="condition"))
    + geom_point(size=1.5)
    + facet_wrap("condition")
    + theme_classic()
)
```

## From AnnData to a publication figure

Publication styling is opt-in and uses the same helper calls. This example uses
Scanpy's bundled PBMC68k subset, so it can be run as written:

```python
import scanpy as sc
import ggann as ag

adata = sc.datasets.pbmc68k_reduced()
genes = ["CD3D", "MS4A1", "NKG7", "GNLY", "CST3"]

with ag.style_context("double-column"):
    embedding = ag.plot_embedding(
        adata, "umap", color="bulk_labels", label=True
    )
    markers = ag.plot_dotplot(
        adata, genes, group_by="bulk_labels", use_raw=True
    )
    figure = ag.compose(
        [embedding, markers],
        ncol=2,
        widths=(0.9, 1.3),
        gap=2,
        guides="keep",
    )

paths = ag.save_publication(
    figure,
    "pbmc_lineages",
    width="double-column",
    height=90,
    formats=("svg", "pdf", "png"),
    dpi=600,
)
```

The context coordinates typography, line weights, palettes, missing-value
colours, embedding axes, and matrix guides. It restores plotnine and Matplotlib
state exactly when the block exits. `save_publication` preserves the requested
canvas instead of tight-cropping it; SVG text stays editable, PDF fonts use the
TrueType path, and raster dimensions follow the requested millimetres and DPI.

The executable [publication workflow](docs/vignettes/publication_panels.md)
extends this to four unequal panels and records the claim, panel map, `n`,
summary definitions, output dimensions, palette, and input fingerprint.

## Installation

Install the public distribution with:

```bash
pip install ggann
```

For development from a checkout:

```bash
python -m pip install -e ".[test,docs]"
```

Optional backends are isolated in extras:

| Extra | Functions |
|---|---|
| `density` | `plot_density` |
| `heatmap` | `plot_clustermap` |
| `upset` | `plot_upset` |
| `pseudobulk` | `pseudobulk` |

```bash
pip install "ggann[density,heatmap,upset,pseudobulk]"
```

## Two interfaces, one data contract

- Use helpers such as `plot_embedding`, `plot_dotplot`, `plot_matrixplot`, and
  `plot_violin` for standard single-cell figures.
- Use `gganndata(...) + ...` when the figure is naturally expressed as a
  plotnine grammar.
- Except for the documented grid-backend functions `plot_clustermap` and
  `plot_upset`, plotting functions return composable plotnine objects.
- Wrap existing helper or bare plotnine construction in `style_context(...)`
  when all panels should share publication defaults. Add later themes, scales,
  coordinates, facets, or annotations normally; the later layer wins.

Bare aesthetic names resolve in this order: `adata.obs`, the selected
expression matrix, then `adata.obsm` coordinates. Use explicit selectors when a
name is ambiguous:

```python
from ggann import aes, gene, gganndata, obs, obsm
from plotnine import geom_point

plot = gganndata(
    adata,
    aes(
        x=obsm("umap", 0),
        y=obsm("umap", 1),
        color=gene("CD3D", use_raw=True),
        shape=obs("condition"),
    ),
) + geom_point()
```

Expression uses `adata.raw` by default when it exists. Set `use_raw=False` for
`.X`, or pass `layer="counts"` for a named layer. A source set on `gene(...)`
overrides the plot-wide source for that gene.

## annplyr interoperability

`ggann` delegates AnnData-to-table extraction to `annplyr` and projects only the
requested expression columns before conversion. Grammar users can set
`gganndata(..., max_matrix_values=...)` to bound cumulative expression and
`obsm` reads. Use `adata.ap.to_df(...)` or `adata.ap.to_tidy(...)` directly when
preparing a custom table with the same annplyr budget.

See the executable [annplyr vignette](docs/vignettes/annplyr_interop.md) and the
[API concepts guide](docs/concepts.md) for source, ownership, sparse, backed,
ordering, and missing-value behavior.

## Performance

Performance claims are accepted only from matched, reproducible comparisons
against Scanpy. The [performance guide](docs/performance.md) defines the timing,
memory, equivalence, and acceptance criteria; it deliberately does not substitute
microbenchmarks for end-to-end results.

## Development checks

```bash
python scripts/run_pyright.py
ruff check src tests benchmarks docs/extensions examples scripts
ruff format --check src tests benchmarks docs/extensions examples scripts
PYTHONPATH=src python -m pytest -q
GGANN_DOCS_OFFLINE=1 sphinx-build -W --keep-going -b html docs docs/_build/html
```

The warning-as-error docs build executes all six scripts in
`examples/vignettes/`. The first-use and publication workflows use Scanpy's
bundled `pbmc68k_reduced` dataset; four focused workflows use a deterministic
PBMC-like fixture. All use a headless renderer, write only to temporary
directories by default, and require no network access.
