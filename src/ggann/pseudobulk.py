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

import warnings

from ._compat import renamed_keyword
from ._grouping import _drop_missing_groups

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
    adata
        Annotated data matrix containing cells and count-like expression.
    sample_col
        obs column identifying the biological replicate (donor / sample / batch).
    group_by
        obs column to split within each sample (e.g. cell type); ``None`` gives one
        profile per sample.
    layer, use_raw
        Where the counts live. decoupler expects **integer counts**; pass the counts
        ``layer=`` (or ``use_raw=True``), or ``skip_checks=True`` to aggregate
        non-counts. ``raw`` is the deprecated spelling of ``use_raw``.
    mode
        Aggregation: ``"sum"`` (default, for count-based DE), ``"mean"`` or
        ``"median"``. Other values are passed to decoupler, which also accepts a
        callable.
    min_cells
        Drop pseudobulk profiles built from fewer than this many cells
        (decoupler records the count in ``obs['psbulk_cells']``). The default of
        10 removes small profiles; a warning reports how many. Pass
        ``min_cells=0`` to keep every profile.
    raw
        Deprecated alias of ``use_raw``.
    skip_checks
        Forwarded to decoupler to permit non-count input.

    Returns
    -------
    A pseudobulk :class:`~anndata.AnnData` whose observations are sample x group
    profiles. It can be passed to ``gganndata`` or expression-summary helpers
    when their required observation columns are present.

    Raises
    ------
    ImportError
        If the ``pseudobulk`` extra is unavailable.
    KeyError
        If ``sample_col``, ``group_by`` or a selected layer is absent.
    ValueError
        If ``layer`` and ``use_raw=True`` are combined, ``use_raw=True`` but
        ``adata.raw`` is absent, ``min_cells`` is negative, or decoupler rejects
        the data.

    Warns
    -----
    UserWarning
        Cells with a missing sample or group label are dropped, and profiles
        below ``min_cells`` are dropped; each warning gives the count.

    Notes
    -----
    A new AnnData object is returned. The input is not mutated; this whole-source
    aggregation is delegated to decoupler and may materialize its selected matrix.
    With ``use_raw=True`` the counts and gene set come from ``adata.raw`` (which
    may hold more genes than ``adata.var``); the result carries the raw genes.

    Examples
    --------
    >>> pb = pseudobulk(adata, sample_col="donor", group_by="cell_type")
    """
    for field in (sample_col, group_by):
        if field is not None and field not in adata.obs:
            raise KeyError(f"{field!r} is not a column of adata.obs.")
    if isinstance(min_cells, bool) or not isinstance(min_cells, int) or min_cells < 0:
        raise ValueError(f"min_cells must be a non-negative integer, got {min_cells!r}.")
    adata = _drop_missing_groups(adata, sample_col, group_by)
    use_raw = renamed_keyword(
        use_raw,
        raw,
        name="use_raw",
        legacy_name="raw",
        default=False,
    )
    if layer is not None and use_raw:
        raise ValueError("Cannot specify use_raw=True and a layer at the same time.")
    if use_raw:
        if adata.raw is None:
            raise ValueError("use_raw=True but adata.raw is None.")
        # decoupler indexes the current .var with raw names, so hand it an
        # AnnData whose matrix and variables are the raw ones (obs kept).
        adata = adata.raw.to_adata()
    dc = _require_decoupler()
    pb = dc.pp.pseudobulk(
        adata,
        sample_col=sample_col,
        groups_col=group_by,
        layer=layer,
        raw=False,
        mode=mode,
        skip_checks=skip_checks,
    )
    if min_cells and "psbulk_cells" in pb.obs.columns:
        keep = (pb.obs["psbulk_cells"] >= min_cells).to_numpy()
        if not keep.all():
            warnings.warn(
                f"dropped {int((~keep).sum())} of {pb.n_obs} pseudobulk profiles built from "
                f"fewer than min_cells={min_cells} cells; pass min_cells=0 to keep them.",
                UserWarning,
                stacklevel=2,
            )
        pb = pb[keep].copy()
    return pb
