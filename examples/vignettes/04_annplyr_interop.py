"""Executable companion to the annplyr interoperability vignette."""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import annplyr as ap
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from _fixture import fingerprint, load_adata
from annplyr.data import X, X_umap, obs, raw
from plotnine import aes, geom_boxplot, geom_col, geom_point, ggplot, labs

import ggann as ag


def main() -> None:
    plain = load_adata(storage="csr")
    before = fingerprint(plain)
    adata = ap.tidy(plain)  # the verbs become methods; ggann accepts it as is

    # ------------------------------------------------------------------ #
    # A tidy pipeline is a ggann input: filter, derive a score, plot it.
    # ------------------------------------------------------------------ #
    scored = adata.filter(obs.percent_mito < 0.03).mutate(
        b_score=(raw.MS4A1 + raw.CD79A) / 2,
    )
    embedding = scored.pipe(ag.plot_embedding, color="b_score")
    embedding.draw()
    grammar = ag.gganndata(scored, aes("UMAP_1", "UMAP_2", color="b_score")) + geom_point()
    grammar.draw()
    # plot_embedding sorts cells so high scores draw on top; align by cell name.
    plotted = embedding.data["b_score"].astype(float).loc[scored.obs_names]
    expected = (
        np.asarray(plain.raw[scored.obs_names, "MS4A1"].X.todense()).ravel()
        + np.asarray(plain.raw[scored.obs_names, "CD79A"].X.todense()).ravel()
    ) / 2
    assert np.allclose(plotted, expected, rtol=1e-5)

    # ------------------------------------------------------------------ #
    # A table from to_wide is a plain data frame: embedding columns are
    # UMAP_1 and UMAP_2, genes keep their bare names.
    # ------------------------------------------------------------------ #
    table = adata.to_wide(
        obs.louvain,
        X_umap,
        X.CST3,
        X.NKG7,
        max_matrix_values=adata.n_obs * 4,  # 2 genes + 2 UMAP columns
    )
    assert set(table.columns) == {"louvain", "UMAP_1", "UMAP_2", "CST3", "NKG7"}
    (
        ggplot(table, aes("UMAP_1", "UMAP_2", color="louvain")) + geom_point() + ag.theme_ggann()
    ).draw()

    # ------------------------------------------------------------------ #
    # join_features(shape="long") keeps .abundance_* as a pandas sparse column.
    # plotnine cannot take box statistics on it; densify the one column.
    # ------------------------------------------------------------------ #
    long = adata.join_features(raw.CST3, raw.NKG7, shape="long")
    long[".abundance_raw"] = long[".abundance_raw"].sparse.to_dense()
    (
        ggplot(long, aes("louvain", ".abundance_raw", fill=".feature"))
        + geom_boxplot(outlier_alpha=0.2)
        + labs(x="", y="log-normalised expression", fill="")
        + ag.theme_publication()
    ).draw()

    # ------------------------------------------------------------------ #
    # Library-size normalisation does not fully remove the depth effect: within
    # the largest myeloid cluster, deeply sequenced cells still show higher
    # apparent CST3. Assert on that cluster only; the rarest ones hold 1-3 cells
    # per depth bin.
    # ------------------------------------------------------------------ #
    summary = adata.summarize(
        mean_CST3=ap.mean(X.CST3),
        by=(obs.louvain, obs.depth),
        max_matrix_values=adata.n_obs,
    )
    if isinstance(summary["mean_CST3"].dtype, pd.SparseDtype):
        summary["mean_CST3"] = summary["mean_CST3"].sparse.to_dense()
    monocytes = summary[summary["louvain"] == "CD14+ Monocytes"].set_index("depth")["mean_CST3"]
    assert monocytes["high"] > monocytes["low"]
    (
        ggplot(summary, aes("louvain", "mean_CST3", fill="depth"))
        + geom_col(position="dodge")
        + labs(x="", y="mean CST3 expression", fill="sequencing depth")
        + ag.theme_publication()
    ).draw()

    # ------------------------------------------------------------------ #
    # count() lists groups in key order; plot_proportions counts the same way.
    # ------------------------------------------------------------------ #
    counts = adata.count(obs.compartment, obs.louvain)
    assert list(counts.columns) == ["compartment", "louvain", "n"]
    proportions = ag.plot_proportions(adata, group_by="louvain", split_by="compartment")
    proportions.draw()
    merged = proportions.data.merge(counts, on=["compartment", "louvain"], suffixes=("", "_ap"))
    assert (merged["n"] == merged["n_ap"]).all()

    # ------------------------------------------------------------------ #
    # A palette stored with annplyr is the palette ggann draws, and filtering
    # keeps colours attached to their labels.
    # ------------------------------------------------------------------ #
    levels = list(adata.obs["louvain"].cat.categories)
    palette = dict(zip(levels, ag.publication_palette("qualitative", len(levels)), strict=True))
    painted = ap.store_palette(adata, "louvain", palette)
    assert ag.obs_colors(painted, "louvain") == palette
    myeloid = painted.filter(obs.compartment == "Myeloid")
    assert ag.obs_colors(myeloid, "louvain") == {
        level: palette[level] for level in myeloid.obs["louvain"].cat.categories
    }
    ag.plot_embedding(myeloid, color="louvain").draw()

    # ggann helpers and direct annplyr use share one object without editing it.
    ag.plot_dotplot(adata, ["CST3", "NKG7"], group_by="louvain", use_raw=False).draw()
    assert fingerprint(plain) == before
    plt.close("all")


if __name__ == "__main__":
    main()
