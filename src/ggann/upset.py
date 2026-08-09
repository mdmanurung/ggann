"""UpSet plots of set intersections, via marsilea.

marsilea's distinctive contribution over ggann's existing heatmaps is the UpSet
plot -- the modern replacement for a Venn diagram when comparing more than a few
sets. In single-cell work the natural sets are marker / DE gene lists per cell
type: an UpSet shows which genes are shared by which combinations of clusters.

Like :func:`ggann.plot_clustermap` (PyComplexHeatmap), this is an escape hatch:
marsilea is an optional dependency and ``plot_upset`` returns marsilea's own
``Upset`` object rather than a plotnine ``ggplot``. Call ``.save(path)`` on it or
let it render to the current matplotlib figure.
"""

from __future__ import annotations

from typing import Iterable, Mapping

from .publication import _active_style

__all__ = ["plot_upset"]


def _require_marsilea():
    try:
        from marsilea.upset import Upset, UpsetData
    except ImportError as exc:  # pragma: no cover - exercised only without the dep
        raise ImportError(
            "plot_upset requires marsilea; install with `pip install marsilea` "
            "(bundled in the ggann[upset] extra)."
        ) from exc
    return Upset, UpsetData


def plot_upset(
    sets: Mapping[str, Iterable[str]],
    *,
    sort_subsets: str = "cardinality",
    min_cardinality: int | None = None,
    min_degree: int | None = None,
    orient: str = "h",
    render: bool = True,
    **kwargs,
):
    """UpSet plot of the intersections between named sets.

    ``sets`` is a mapping of ``{set_name: iterable_of_members}`` -- e.g. the top
    marker genes of each cell type (build it from :func:`ggann.rank_genes_df`) to
    see which markers are shared across clusters.

    ``sort_subsets`` orders the intersection bars (``"cardinality"`` or
    ``"degree"``); ``min_cardinality`` / ``min_degree`` drop small or
    low-order intersections. Extra ``kwargs`` pass through to marsilea's ``Upset``.
    Returns the marsilea ``Upset`` (rendered to the current figure when
    ``render=True``).

    Parameters
    ----------
    sets : mapping of str to iterable
        Named member sets.
    sort_subsets : {"cardinality", "degree"}
        Intersection ordering.
    min_cardinality : int, optional
        Minimum intersection size.
    min_degree : int, optional
        Minimum number of participating sets.
    orient : {"h", "v"}
        Backend orientation.
    render : bool
        Render immediately.
    **kwargs
        Passed to ``marsilea.Upset``.

    Returns
    -------
    marsilea.Upset
        Grid-backend object; this is not a plotnine object.

    Raises
    ------
    ImportError
        If the ``upset`` extra is unavailable.
    TypeError
        If ``sets`` is not a mapping.
    ValueError
        If fewer than two sets are supplied or filters are invalid.

    Notes
    -----
    Input iterables are copied into membership data; no AnnData is involved.

    Examples
    --------
    >>> up = plot_upset({"T": {"CD3D", "IL7R"}, "NK": {"NKG7", "GNLY"}})
    """
    Upset, UpsetData = _require_marsilea()

    if not isinstance(sets, Mapping):
        raise TypeError("sets must be a mapping of {set_name: iterable_of_members}.")
    names = list(sets.keys())
    if len(names) < 2:
        raise ValueError("plot_upset needs at least two sets to intersect.")
    members = [set(sets[name]) for name in names]

    data = UpsetData.from_sets(members, sets_names=names)
    upset = Upset(
        data,
        orient=orient,
        sort_subsets=sort_subsets,
        min_cardinality=min_cardinality,
        min_degree=min_degree,
        **kwargs,
    )
    if render:
        upset.render()
    if style := _active_style():
        upset._ggann_publication_style = style
    return upset
