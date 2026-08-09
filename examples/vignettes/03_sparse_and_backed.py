"""Executable companion to the sparse and backed AnnData vignette."""

from __future__ import annotations

import tempfile
from pathlib import Path

import matplotlib
from anndata import read_h5ad

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from _fixture import fingerprint, make_adata

import ggann as ag


def main() -> None:
    for storage in ("csr", "csc"):
        adata = make_adata(storage=storage)
        before = fingerprint(adata)
        ag.plot_dotplot(
            adata,
            ["CD3D", "NKG7", "MS4A1", "CST3"],
            group_by="cell_type",
            use_raw=False,
        ).draw()
        assert fingerprint(adata) == before
        plt.close("all")

    with tempfile.TemporaryDirectory(prefix="ggann-vignette-") as tmp:
        for storage in ("dense", "csr", "csc"):
            path = Path(tmp) / f"fixture-{storage}.h5ad"
            make_adata(storage=storage).write_h5ad(path)
            backed = read_h5ad(path, backed="r")
            try:
                assert backed.isbacked
                ag.plot_embedding(backed, "umap", color="cell_type").draw()
                ag.plot_matrixplot(
                    backed,
                    ["CD3D", "NKG7", "MS4A1", "CST3"],
                    group_by="cell_type",
                    use_raw=False,
                ).draw()
            finally:
                backed.file.close()
            plt.close("all")


if __name__ == "__main__":
    main()
