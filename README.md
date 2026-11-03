# ggann

`ggann` turns `AnnData` into composable plotnine graphics and concise
single-cell plotting helpers. Start with a one-call biological summary, refine
it with the grammar of graphics, and export an editable multi-panel figure at
its final physical size.

[Documentation](https://mdmanurung.github.io/ggann/) ·
[Coming from Scanpy](https://mdmanurung.github.io/ggann/from_scanpy.html) ·
[Custom figures](https://mdmanurung.github.io/ggann/grammar.html) ·
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
the 2,638-cell PBMC3K dataset from Scanpy's and Seurat's clustering tutorials,
so it can be run as written:

```python
import scanpy as sc
import ggann as ag

adata = sc.datasets.pbmc3k_processed()  # downloads once, then caches
genes = ["CD3D", "MS4A1", "NKG7", "GNLY", "CST3"]

with ag.style_context("double-column"):
    embedding = ag.plot_embedding(
        adata, "umap", color="louvain", label=True
    )
    markers = ag.plot_dotplot(
        adata, genes, group_by="louvain", use_raw=True
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
TrueType path, and raster dimensions follow the requested millimetres and DPI. For
`plot_clustermap` and `tidy_heatmap` figures, it also shrinks the heatmap body
until every label and legend fits the canvas.

The [publication guide](docs/publication.md) extends this to a four-panel figure.

## Installation

ggann needs `annplyr>=0.5,<0.6` and `plotnine-extra>=0.3.1,<0.4`. Neither release
is on PyPI yet, so `pip install ggann` does not resolve today. The package
metadata cannot name Git URLs, so install both from GitHub first, at the commits
CI tests against:

```bash
pip install "annplyr @ git+https://github.com/mdmanurung/annplyr@ae2db5402c8b1c714bb12906ff82a40fbb0d8de2"
pip install "plotnine-extra @ git+https://github.com/mdmanurung/plotnine-extra@fcbd88f64760a2dd6e8200028a1a00909aa2c9f8"
pip install "ggann @ git+https://github.com/mdmanurung/ggann"
```

For development from a checkout, install the two dependencies as above, then:

```bash
python -m pip install -e ".[test,docs]"
```

Once ggann 0.1.0 is released and both dependencies are on PyPI, this will work:

```bash
pip install ggann
```

Optional backends are isolated in extras:

| Extra | Functions |
|---|---|
| `density` | `plot_density` |
| `heatmap` | `plot_clustermap`, `tidy_heatmap` |
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
- Except for the documented grid-backend functions `plot_clustermap`,
  `tidy_heatmap` and `plot_upset`, plotting functions return composable
  plotnine objects.
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

`ggann` depends on `annplyr` 0.5 and reads `AnnData` through it. An annplyr
result plots without an adapter: `ap.tidy(adata)` returns a `TidyAnnData`
(an `AnnData` subclass), and `gganndata` and every `plot_*` helper accept it.

```python
import annplyr as ap
from annplyr.data import X_umap, obs, raw
from plotnine import aes, geom_point, ggplot

import ggann as ag

adata = ap.tidy(adata)
scored = adata.filter(obs.percent_mito < 0.05).mutate(b_score=(raw.MS4A1 + raw.CD79A) / 2)

scored.pipe(ag.plot_embedding, color="b_score")  # helper, chained
ag.gganndata(scored, aes("UMAP_1", "UMAP_2", color="b_score")) + geom_point()  # grammar

table = adata.to_wide(obs.louvain, X_umap)  # columns: louvain, UMAP_1, UMAP_2
ggplot(table, aes("UMAP_1", "UMAP_2", color="louvain")) + geom_point() + ag.theme_ggann()
```

Tables from `to_wide`, `join_features`, `summarize` and `count` are DataFrames that
plotnine draws directly. `ap.store_palette` colours are the ones `plot_embedding` and
`scale_color_obs` draw.

`ggann` projects only the requested expression columns before conversion. Use
`gganndata(..., max_matrix_values=...)` or the same keyword on an annplyr verb to
bound what is read.

See the [custom figures guide](docs/grammar.md) for annplyr pipelines, name
resolution and large-data budgets.

## Performance

Performance claims are accepted only from matched, reproducible comparisons
against Scanpy. The [notes page](docs/notes.md) gives the timing,
memory, equivalence, and acceptance criteria; it deliberately does not substitute
microbenchmarks for end-to-end results.

## Citing

See [Citing ggann](docs/notes.md#citing) for how to cite ggann and the tools behind each figure.

## Development checks

```bash
python scripts/run_pyright.py
ruff check src tests benchmarks docs/extensions examples scripts
ruff format --check src tests benchmarks docs/extensions examples scripts
PYTHONPATH=src python -m pytest -q
GGANN_DOCS_OFFLINE=1 sphinx-build -W --keep-going -b html docs docs/_build/html
```

The warning-as-error docs build executes all six scripts in
`examples/vignettes/`. All six run on the same real dataset,
`scanpy.datasets.pbmc3k_processed`. Fetch it once with `python
scripts/fetch_datasets.py`; the scripts then use a headless renderer, write
only to temporary directories by default, and reach no network themselves.
