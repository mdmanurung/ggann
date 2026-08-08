"""Pseudobulk aggregation, so cell-level AnnData can be summarised to sample level.

Aggregates single cells to one profile per ``sample_col`` x ``group_by`` (e.g. per
donor per cell type) with decoupler's ``pp.pseudobulk`` -- the same aggregation
liana's MOFA-cellular workflow builds on. It returns a **new AnnData** that can
be plotted when it retains the observation columns required by the chosen call:

    pb = ag.pseudobulk(adata, sample_col="donor", group_by="cell_type")
    ag.plot_dotplot(pb, markers, "cell_type", use_raw=False)
    gganndata(pb, aes("cell_type", gene("CD3D"), fill="condition")) + geom_boxplot()

decoupler is an optional dependency (``ggann[pseudobulk]``).
"""

from __future__ import annotations

from ._compat import renamed_keyword

__all__ = ["pseudobulk"]


def _require_decoupler():
    try:
        import decoupler as dc
    except ImportError as exc:  # pragma: no cover - only without the dep
        raise ImportError(
            "pseudobulk requires decoupler; install with `pip install decoupler` "
            "(bundled in the ggann[pseudobulk] extra)."
        ) from exc
    return dc


def pseudobulk(
    adata,
    sample_col: str,
    group_by: str | None = None,
    *,
    layer: str | None = None,
    use_raw: bool | None = None,
    raw: bool | None = None,
    mode: str = "sum",
    min_cells: int = 10,
    skip_checks: bool = False,
):
    """Aggregate cells to a pseudobulk AnnData (one profile per sample x group).

    Parameters
    ----------
    sample_col
        obs column identifying the biological replicate (donor / sample / batch).
    group_by
        obs column to split within each sample (e.g. cell type); ``None`` gives one
        profile per sample.
    layer / use_raw
        Where the counts live. decoupler expects **integer counts**; pass the counts
        ``layer=`` (or ``use_raw=True``), or ``skip_checks=True`` to aggregate
        non-counts. ``raw`` is the deprecated spelling of ``use_raw``.
    mode
        Aggregation: ``"sum"`` (default, for count-based DE) or ``"mean"``.
    min_cells
        Drop pseudobulk profiles built from fewer than this many cells
        (decoupler records the count in ``obs['psbulk_cells']``).

    Returns
    -------
    A pseudobulk :class:`~anndata.AnnData` whose observations are sample x group
    profiles. It can be passed to ``gganndata`` or expression-summary helpers
    when their required observation columns are present.
    """
    use_raw = renamed_keyword(
        use_raw,
        raw,
        name="use_raw",
        legacy_name="raw",
        default=False,
    )
    if layer is not None and use_raw:
        raise ValueError("Cannot specify use_raw=True and a layer at the same time.")
    dc = _require_decoupler()
    pb = dc.pp.pseudobulk(
        adata,
        sample_col=sample_col,
        groups_col=group_by,
        layer=layer,
        raw=use_raw,
        mode=mode,
        skip_checks=skip_checks,
    )
    if min_cells and "psbulk_cells" in pb.obs.columns:
        pb = pb[pb.obs["psbulk_cells"] >= min_cells].copy()
    return pb
