"""Executable companion to the annplyr interoperability vignette."""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import annplyr as ap
import matplotlib.pyplot as plt
import pandas as pd
from _fixture import fingerprint, load_adata
from plotnine import aes, geom_col, ggplot, labs

import ggann as ag


def main() -> None:
    adata = load_adata(storage="csr")
    before = fingerprint(adata)

    # annplyr projects just two genes and rejects a request above this budget.
    table = adata.ap.to_df(
        obs=["louvain", "depth"],
        x=["CST3", "NKG7"],
        max_matrix_values=2 * adata.n_obs,
    )
    assert {"louvain", "depth", "CST3", "NKG7"} <= set(table.columns)

    summary = adata.ap.summarize(
        x={"mean_CST3": ap.mean(ap.col("CST3"))},
        by=["louvain", "depth"],
        max_matrix_values=adata.n_obs,
    )
    if isinstance(summary["mean_CST3"].dtype, pd.SparseDtype):
        summary["mean_CST3"] = summary["mean_CST3"].sparse.to_dense()

    # Library-size normalisation does not fully remove the depth effect: within
    # the largest myeloid cluster, deeply sequenced cells still show higher
    # apparent CST3. Assert on that cluster only; the rarest ones hold 1-3 cells
    # per depth bin.
    monocytes = summary[summary["louvain"] == "CD14+ Monocytes"].set_index("depth")["mean_CST3"]
    assert monocytes["high"] > monocytes["low"]

    custom = (
        ggplot(summary, aes("louvain", "mean_CST3", fill="depth"))
        + geom_col(position="dodge")
        + labs(x="", y="mean CST3 expression", fill="sequencing depth")
        + ag.theme_publication()
    )
    custom.draw()
    ag.plot_dotplot(
        adata,
        ["CST3", "NKG7"],
        group_by="louvain",
        use_raw=False,
    ).draw()
    assert fingerprint(adata) == before
    plt.close("all")


if __name__ == "__main__":
    main()
