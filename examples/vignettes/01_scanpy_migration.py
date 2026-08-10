"""Build a reusable PBMC marker figure from a familiar Scanpy workflow."""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import scanpy as sc
from _fixture import fingerprint, load_adata
from plotnine import ggplot
from plotnine.composition import Compose

import ggann as ag


def main() -> None:
    adata = load_adata()
    # CD3D and LYZ did not survive this dataset's highly-variable-gene
    # selection, so the marker panel reads them from ``.raw`` with use_raw=True.
    genes = ["CD3D", "IL7R", "MS4A1", "NKG7", "GNLY", "LYZ"]

    # A familiar Scanpy call is a useful starting point for the translation.
    # Use a copy because third-party plotting is outside ggann's ownership
    # contract and may cache plotting metadata on AnnData.
    sc.pl.embedding(adata.copy(), basis="umap", color="louvain", show=False)
    plt.close("all")

    before = fingerprint(adata)
    embedding = ag.plot_embedding(
        adata,
        basis="umap",
        color="louvain",
        label=True,
    )
    markers = ag.plot_dotplot(
        adata,
        genes,
        group_by="louvain",
        use_raw=True,
    )
    distribution = ag.plot_violin(
        adata,
        ["CD3D"],
        group_by="louvain",
        use_raw=True,
        add_box=True,
        stats=False,
    )

    figure = ag.compose(
        [embedding, markers],
        ncol=2,
        widths=(0.9, 1.3),
        gap=2,
        tag_levels=None,
    )

    assert all(isinstance(plot, ggplot) for plot in (embedding, markers, distribution))
    assert isinstance(figure, Compose)
    assert len(embedding.data) == adata.n_obs
    assert len(distribution.data) == adata.n_obs
    assert set(markers.data["feature"]) == set(genes)
    figure.draw(show=False)
    distribution.draw(show=False)
    assert fingerprint(adata) == before
    plt.close("all")


if __name__ == "__main__":
    main()
