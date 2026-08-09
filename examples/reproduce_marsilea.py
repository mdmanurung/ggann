"""Generate ggann versions of the scanpy and marsilea plotting tutorial.

Run ``python examples/reproduce_marsilea.py`` with the heatmap extra installed.
Images are written to ``examples/_output``.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import scanpy as sc

import ggann as ag

OUT = Path(__file__).parent / "_output"


def _save(plot, name: str, **kwargs) -> None:
    path = OUT / name
    plot.save(path, verbose=False, **kwargs)
    print("wrote", path)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    adata = sc.datasets.pbmc68k_reduced()
    group = "bulk_labels"
    markers = ["CD3D", "CD8A", "NKG7", "GNLY", "MS4A1", "FCGR3A", "CST3"]
    markers = [gene for gene in markers if gene in adata.raw.var_names]

    _save(
        ag.plot_embedding(adata, basis="umap", color=group),
        "01_umap_clusters.png",
        width=7,
        height=5,
        dpi=120,
    )
    _save(
        ag.plot_embedding(adata, basis="umap", color="CD3D"),
        "02_umap_CD3D.png",
        width=6,
        height=5,
        dpi=120,
    )
    _save(
        ag.plot_dotplot(adata, markers, group),
        "03_dotplot.png",
        width=7,
        height=5,
        dpi=120,
    )
    _save(
        ag.plot_matrixplot(adata, markers, group, standard_scale="var"),
        "04_matrixplot.png",
        width=7,
        height=5,
        dpi=120,
    )
    _save(
        ag.plot_violin(adata, markers, group),
        "05_violin.png",
        width=6,
        height=12,
        dpi=120,
    )

    plt.figure(figsize=(7, 5))
    ag.plot_clustermap(
        adata,
        markers,
        group_by=group,
        standard_scale="var",
        z_score=None,
    )
    clustermap_path = OUT / "06_clustermap.png"
    plt.savefig(clustermap_path, dpi=120, bbox_inches="tight")
    plt.close("all")
    print("wrote", clustermap_path)


if __name__ == "__main__":
    main()
