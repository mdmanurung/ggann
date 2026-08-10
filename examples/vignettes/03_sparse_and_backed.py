"""Executable companion to the sparse and backed AnnData vignette."""

from __future__ import annotations

import tempfile
from pathlib import Path

import matplotlib
from anndata import read_h5ad

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from _fixture import MARKERS, fingerprint, load_adata

import ggann as ag


def main() -> None:
    for storage in ("csr", "csc"):
        adata = load_adata(storage=storage)
        before = fingerprint(adata)
        ag.plot_dotplot(
            adata,
            MARKERS,
            group_by="louvain",
            use_raw=False,
        ).draw()
        assert fingerprint(adata) == before
        plt.close("all")

    with tempfile.TemporaryDirectory(prefix="ggann-vignette-") as tmp:
        for storage in ("dense", "csr", "csc"):
            path = Path(tmp) / f"pbmc3k-{storage}.h5ad"
            load_adata(storage=storage).write_h5ad(path)
            backed = read_h5ad(path, backed="r")
            try:
                assert backed.isbacked
                ag.plot_embedding(backed, "umap", color="louvain").draw()
                ag.plot_matrixplot(
                    backed,
                    MARKERS,
                    group_by="louvain",
                    use_raw=False,
                ).draw()
            finally:
                backed.file.close()
            plt.close("all")


if __name__ == "__main__":
    main()
