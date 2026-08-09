"""Executable companion to the publication multi-panel vignette."""

from __future__ import annotations

import tempfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

from _fixture import make_adata

import ggann as ag


def main() -> None:
    adata = make_adata()
    genes = ["CD3D", "NKG7", "MS4A1", "CST3"]
    panels = [
        ag.plot_embedding(adata, "umap", color="cell_type", label=True),
        ag.plot_dotplot(adata, genes, group_by="cell_type"),
        ag.plot_violin(adata, ["CD3D"], group_by="cell_type"),
        ag.plot_proportions(adata, "cell_type", split_by="condition", normalize=True),
    ]
    figure = ag.compose(panels, ncol=2, tag_levels="A")

    with tempfile.TemporaryDirectory(prefix="ggann-vignette-") as tmp:
        output = Path(tmp) / "figure.png"
        figure.save(output, width=180, height=140, units="mm", dpi=90, verbose=False)
        assert output.stat().st_size > 0


if __name__ == "__main__":
    main()
