"""From exploratory AnnData plots to a publication-ready PBMC figure.

The scientific claim is that broad PBMC lineages occupy distinct UMAP
neighbourhoods and show coherent lineage-marker programs. The four panels map
that claim to position, marker expression, a log-normalized expression
distribution, and composition across sequencing-depth halves. Run with
``--output figures/pbmc``
to retain the figure and its machine-readable manifest; the documentation
build uses a temporary directory.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import platform
import struct
import tempfile
from contextlib import nullcontext
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import pandas as pd
from _fixture import fingerprint, load_adata
from plotnine import labs, theme

import ggann as ag
from ggann._palette_qa import palette_accessibility_report

CLAIM = (
    "Broad PBMC lineages occupy distinct UMAP neighbourhoods and show "
    "coherent lineage-marker programs."
)
GENES = ["CD3D", "MS4A1", "NKG7", "GNLY", "CST3"]
LINEAGES = ["T cell", "B cell", "NK cell", "Monocyte", "Dendritic", "Megakaryocyte"]
PANEL_MAP = {
    "a": "all-cell UMAP, coloured and directly labelled by broad lineage",
    "b": "mean log-normalized raw expression and fraction detected for five lineage markers",
    "c": "all-cell CD3D distributions with median and interquartile range",
    "d": "lineage proportions within each half of the sequencing-depth distribution",
}


def load_data():
    """Load the processed PBMC3K dataset and add a broad-lineage annotation."""
    adata = load_adata()
    lineage = (
        adata.obs["louvain"]
        .astype("string")
        .map(
            {
                "CD4 T cells": "T cell",
                "CD8 T cells": "T cell",
                "B cells": "B cell",
                "NK cells": "NK cell",
                "CD14+ Monocytes": "Monocyte",
                "FCGR3A+ Monocytes": "Monocyte",
                "Dendritic cells": "Dendritic",
                "Megakaryocytes": "Megakaryocyte",
            }
        )
    )
    if lineage.isna().any():
        raise AssertionError("The PBMC3K Louvain labels changed; review the lineage mapping.")
    adata.obs["lineage"] = pd.Categorical(lineage, categories=LINEAGES, ordered=True)

    # One explicit vocabulary is reused by every panel that represents lineage.
    palette = ag.publication_palette("qualitative", categories=LINEAGES)
    adata.uns["lineage_colors"] = [palette[lineage] for lineage in LINEAGES]
    return adata, palette


def build_panels(adata):
    """Build one panel per part of the claim using the same ergonomic helper API."""
    return [
        ag.plot_embedding(
            adata,
            "umap",
            color="lineage",
            label=True,
            size=1.0,
        )
        + labs(title="PBMC lineage map", color="lineage"),
        ag.plot_dotplot(
            adata,
            GENES,
            group_by="lineage",
            use_raw=True,
        )
        + labs(title="Lineage-marker programs"),
        ag.plot_violin(
            adata,
            ["CD3D"],
            group_by="lineage",
            use_raw=True,
            add_box=True,
            add_points=False,
            stats=False,
        )
        + labs(title="CD3D distributions")
        + theme(legend_position="none"),
        ag.plot_proportions(
            adata,
            group_by="lineage",
            split_by="depth",
            normalize=True,
        )
        + labs(title="Composition by sequencing depth"),
    ]


def _assert_same_prepared_content(before, after) -> None:
    """Prove that publication styling did not change observations or summaries."""
    for exploratory, publication in zip(before, after, strict=True):
        pd.testing.assert_frame_equal(
            exploratory.data,
            publication.data,
            check_exact=True,
            check_categorical=True,
        )


def _png_size(path: Path) -> tuple[int, int]:
    """Read the PNG IHDR dimensions without adding an image-library dependency."""
    header = path.read_bytes()[:24]
    if header[:8] != b"\x89PNG\r\n\x1a\n":
        raise AssertionError(f"{path} is not a PNG file.")
    return struct.unpack(">II", header[16:24])


def _sha256(path: Path) -> str:
    """Hash a retained output without loading it into a plotting library."""
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def render(output: Path) -> dict[str, object]:
    """Render before/after figures and return the evidence manifest."""
    output.mkdir(parents=True, exist_ok=True)
    adata, palette = load_data()
    input_sha256 = fingerprint(adata)

    exploratory_panels = build_panels(adata)
    exploratory_figure = ag.compose(
        exploratory_panels,
        ncol=2,
        widths=(0.95, 1.25),
        heights=(1.0, 1.0),
        gap=2,
        tag_levels="auto",
    )
    before_paths = ag.save_publication(
        exploratory_figure,
        output / "pbmc_exploratory.png",
        width="double-column",
        height=120,
        dpi=300,
    )

    with ag.style_context("double-column", dpi=300) as style:
        publication_panels = build_panels(adata)
        publication_figure = ag.compose(
            publication_panels,
            ncol=2,
            widths=(0.95, 1.25),
            heights=(1.0, 1.0),
            gap=2,
            guides="keep",
            tag_levels="auto",
        ) + theme(legend_position="right")

    _assert_same_prepared_content(exploratory_panels, publication_panels)
    final_paths = ag.save_publication(
        publication_figure,
        output / "pbmc_publication",
        width="double-column",
        height=120,
        formats=("svg", "pdf", "png", "tiff"),
        dpi=300,
        background="white",
    )

    expected_pixels = (round(183 / 25.4 * 300), round(120 / 25.4 * 300))
    publication_png = output / "pbmc_publication.png"
    if _png_size(publication_png) != expected_pixels:
        raise AssertionError("The publication PNG does not match 183 x 120 mm at 300 DPI.")
    svg = (output / "pbmc_publication.svg").read_text()
    if "<text" not in svg or "<image" in svg:
        raise AssertionError("Default SVG output must contain vector text and no image layer.")
    if fingerprint(adata) != input_sha256:
        raise AssertionError("Plot construction or export mutated the AnnData input.")

    accessibility = palette_accessibility_report(tuple(palette.values()))
    separations = accessibility["minimum_ciede2000"]
    accessibility["thresholds"] = {"normal": 10.0, "colour_vision_simulations": 5.0}
    accessibility["passes"] = bool(
        separations["normal"] >= 10
        and all(separations[name] >= 5 for name in ("protanopia", "deuteranopia", "tritanopia"))
    )
    (output / "accessibility.json").write_text(
        json.dumps(accessibility, indent=2, sort_keys=True) + "\n"
    )

    retained_paths = [*before_paths, *final_paths, output / "accessibility.json"]

    n_by_lineage = adata.obs["lineage"].value_counts(sort=False)
    manifest: dict[str, object] = {
        "dataset": "scanpy.datasets.pbmc3k_processed",
        "claim": CLAIM,
        "panel_map": PANEL_MAP,
        "n_cells": int(adata.n_obs),
        "n_by_lineage": {str(key): int(value) for key, value in n_by_lineage.items()},
        "statistics": {
            "a": {"n": int(adata.n_obs), "center": None, "interval": None, "test": None},
            "b": {
                "n": int(adata.n_obs),
                "center": "arithmetic mean of log-normalized raw expression",
                "interval": None,
                "test": None,
                "size": "fraction of cells with expression above zero",
            },
            "c": {
                "n": int(adata.n_obs),
                "center": "median",
                "interval": "interquartile range",
                "test": "none; descriptive panel",
            },
            "d": {
                "n": int(adata.n_obs),
                "center": None,
                "interval": None,
                "test": "none; descriptive panel",
                "denominator": "all cells within each half of the sequencing-depth distribution",
            },
        },
        "palette": palette,
        "palette_accessibility": accessibility,
        "style": style.to_dict(),
        "input_sha256": input_sha256,
        "canvas_mm": [183, 120],
        "dpi": 300,
        "raster_pixels": list(expected_pixels),
        "exploratory_outputs": [path.name for path in before_paths],
        "publication_outputs": [path.name for path in final_paths],
        "output_sha256": {path.name: _sha256(path) for path in retained_paths},
        "software": {
            "python": platform.python_version(),
            "ggann": ag.__version__,
            "anndata": importlib.metadata.version("anndata"),
            "matplotlib": importlib.metadata.version("matplotlib"),
            "pandas": importlib.metadata.version("pandas"),
            "plotnine": importlib.metadata.version("plotnine"),
            "scanpy": importlib.metadata.version("scanpy"),
        },
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        help="Directory in which to retain SVG, PDF, PNG, TIFF, and manifest outputs.",
    )
    args = parser.parse_args()
    temporary = (
        tempfile.TemporaryDirectory(prefix="ggann-publication-vignette-")
        if args.output is None
        else nullcontext(str(args.output))
    )
    with temporary as directory:
        manifest = render(Path(directory))
        assert manifest["n_cells"] == 2638


if __name__ == "__main__":
    main()
