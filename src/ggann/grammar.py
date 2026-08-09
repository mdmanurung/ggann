"""The plotnine-native entrypoint: ``gganndata(adata) + aes(...) + geom_*()``.

``gganndata`` eagerly resolves the names referenced in the aesthetic mapping into
a tidy DataFrame (via annplyr) and returns a *real* :class:`plotnine.ggplot`. It
deliberately does **not** subclass ``ggplot`` or intercept plotnine's draw-time
resolution. Because it returns a plain ggplot, existing plotnine geoms, stats,
scales, facets, and themes compose without adapters.
"""

from __future__ import annotations

from plotnine import aes, ggplot

from ._resolve import ObsmRef, Ref, plain_name, resolve_frame
from .publication import _family_theme

__all__ = ["gganndata", "aes"]


def _referenced_names(mapping) -> list:
    """Extract the values an ``aes()`` mapping refers to, for resolution.

    ``aes`` is a plain dict subclass, so its values are the mapped expressions.
    We keep strings (bare or prefixed like ``"gene:CD3D"``) and
    :class:`Ref`/:class:`ObsmRef` accessors; anything else (stages,
    ``after_stat``, callables) is assumed to reference already-present columns.
    """
    if mapping is None:
        return []
    return [v for v in mapping.values() if isinstance(v, (str, Ref, ObsmRef))]


def _plain_mapping(adata, mapping):
    """Rewrite accessor / prefix-string aes values to their plain column names.

    plotnine evaluates aes values against the DataFrame, so ``"gene:CD3D"``,
    ``obsm("umap", 0)`` etc. must be swapped for the string they resolve to
    (``"CD3D"``, ``"UMAP_1"``) before the mapping reaches ``ggplot``. The
    DataFrame from :func:`resolve_frame` is keyed by exactly those names.
    """
    plain = aes()
    for key, value in mapping.items():
        if isinstance(value, (Ref, ObsmRef)):
            plain[key] = plain_name(adata, value)
        elif isinstance(value, str):
            plain[key] = plain_name(adata, value)
        else:
            plain[key] = value
    return plain


def gganndata(
    adata,
    mapping=None,
    *,
    layer: str | None = None,
    use_raw: bool | None = None,
    max_matrix_values: int | None = None,
    add_theme: bool = True,
):
    """Start a ggplot over an ``AnnData``.

    Parameters
    ----------
    adata
        The annotated data matrix.
    mapping
        A plotnine ``aes(...)``. Names are resolved across ``obs`` columns, genes
        (in ``X`` / a ``layer`` / ``raw``) and embedding coordinates such as
        ``"UMAP_1"``. Use :func:`ggann.gene` / :func:`ggann.obs` to force a source.
    layer
        Read expression from ``adata.layers[layer]`` instead of ``X``.
    use_raw
        Read expression from ``adata.raw``. Defaults to ``True`` when ``adata.raw``
        exists and no ``layer`` is given (scanpy's convention).
    max_matrix_values
        Maximum cumulative number of expression and ``obsm`` values that may be
        materialized while resolving the mapping. The request is rejected before
        extraction when it exceeds this annplyr-compatible budget. Observation
        metadata does not count toward the budget. ``None`` disables the limit.
    add_theme
        Add :func:`theme_ggann` to the returned plot.

    Returns
    -------
    plotnine.ggplot
        Compose it with any plotnine layer, e.g.
        ``gganndata(adata, aes("UMAP_1", "UMAP_2", color="louvain")) + geom_point()``.

    Raises
    ------
    KeyError
        If an explicitly selected observation, gene, embedding, or layer is absent.
    ValueError
        If ``layer`` and ``use_raw=True`` are combined.
    annplyr.AnnplyrError
        If ``max_matrix_values`` is negative or the cumulative projection exceeds it.

    Notes
    -----
    Only mapped columns are projected through annplyr. The input ``adata`` is not
    mutated, and the returned object is an ordinary composable plotnine plot.

    Examples
    --------
    >>> from plotnine import geom_point
    >>> p = gganndata(adata, aes("UMAP_1", "UMAP_2", color="cell_type"))
    >>> p = p + geom_point()
    """
    mapping = aes() if mapping is None else mapping
    df = resolve_frame(
        adata,
        _referenced_names(mapping),
        layer=layer,
        use_raw=use_raw,
        max_matrix_values=max_matrix_values,
    )
    plot = ggplot(df, _plain_mapping(adata, mapping))
    if add_theme:
        plot = plot + _family_theme("standard")
    return plot
