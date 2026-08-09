"""Build a reusable PBMC marker figure from a familiar Scanpy workflow."""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import scanpy as sc
from _fixture import fingerprint
from plotnine import ggplot
from plotnine.composition import Compose

import ggann as ag


def main() -> None:
    adata = sc.datasets.pbmc68k_reduced()
    genes = ["CD3D", "MS4A1", "NKG7", "GNLY", "CST3"]

    # A familiar Scanpy call is a useful starting point for the translation.
    # Use a copy because third-party plotting is outside ggann's ownership
    # contract and may cache plotting metadata on AnnData.
    sc.pl.embedding(adata.copy(), basis="umap", color="bulk_labels", show=False)
    plt.close("all")

    before = fingerprint(adata)
    embedding = ag.plot_embedding(
        adata,
        basis="umap",
        color="bulk_labels",
        label=True,
    )
    markers = ag.plot_dotplot(
        adata,
        genes,
        group_by="bulk_labels",
        use_raw=True,
    )
    distribution = ag.plot_violin(
        adata,
        ["CD3D"],
        group_by="bulk_labels",
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
