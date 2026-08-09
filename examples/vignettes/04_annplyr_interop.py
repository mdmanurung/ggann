"""Executable companion to the annplyr interoperability vignette."""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import annplyr as ap
import matplotlib.pyplot as plt
import pandas as pd
from _fixture import fingerprint, make_adata
from plotnine import aes, geom_col, ggplot, labs

import ggann as ag


def main() -> None:
    adata = make_adata(storage="csr")
    before = fingerprint(adata)

    # annplyr projects just two genes and rejects a request above this budget.
    table = adata.ap.to_df(
        obs=["cell_type", "condition"],
        x=["CD3D", "NKG7"],
        max_matrix_values=2 * adata.n_obs,
    )
    assert {"cell_type", "condition", "CD3D", "NKG7"} <= set(table.columns)

    summary = adata.ap.summarize(
        x={"mean_CD3D": ap.mean(ap.col("CD3D"))},
        by=["cell_type", "condition"],
        max_matrix_values=adata.n_obs,
    )
    if isinstance(summary["mean_CD3D"].dtype, pd.SparseDtype):
        summary["mean_CD3D"] = summary["mean_CD3D"].sparse.to_dense()
    custom = (
        ggplot(summary, aes("cell_type", "mean_CD3D", fill="condition"))
        + geom_col(position="dodge")
        + labs(x="", y="mean CD3D expression", fill="condition")
        + ag.theme_publication()
    )
    custom.draw()
    ag.plot_dotplot(
        adata,
        ["CD3D", "NKG7"],
        group_by="cell_type",
        use_raw=False,
    ).draw()
    assert fingerprint(adata) == before
    plt.close("all")


if __name__ == "__main__":
    main()
