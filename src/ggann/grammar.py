"""The plotnine-native entrypoint: ``gganndata(adata) + aes(...) + geom_*()``.

``gganndata`` resolves the AnnData fields named in the aesthetic mapping into a
tidy DataFrame (via annplyr) and returns a *real* :class:`plotnine.ggplot`
subclass. Fields named *later* (``+ aes(...)``, ``+ geom_*(aes(...))``,
``+ facet_*(...)``) are read from the same AnnData when they are added, within
the same materialization budget. Because it is a ggplot, existing plotnine
geoms, stats, scales, facets, and themes compose without adapters.
"""

from __future__ import annotations

import ast
from copy import copy
from dataclasses import dataclass, replace

import pandas as pd
from patsy.builtins import Q
from plotnine import aes, ggplot
from plotnine.mapping._env import Environment

from ._resolve import ObsmRef, Ref, field_names, parse_token, plain_name, resolve_frame
from .publication import _family_theme

__all__ = ["gganndata", "aes"]


# --------------------------------------------------------------------------- #
# Which AnnData fields does an aes value need?
# --------------------------------------------------------------------------- #
def _identifiers(expression: str) -> list[str]:
    """Field identifiers of an aes expression (the supported subset).

    Supported: a bare name, or a Python expression over names, numbers, strings,
    operators, calls and attribute access, e.g. ``"np.log1p(n_counts) / 2"`` or
    ``'factor(louvain)'``. A name is a *field* unless it is called
    (``factor(...)``) or is the base of an attribute (``np.log1p``). A gene that
    is not a Python identifier is written ``Q("HLA-DRA")`` or passed whole as
    ``gene("HLA-DRA")``. Nothing is evaluated here; plotnine evaluates the
    expression against the plot data, in the caller's namespace.
    """
    try:
        tree = ast.parse(expression.strip(), mode="eval")
    except SyntaxError:
        return []
    skip: set[int] = set()
    found: list[str] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Call):
            if isinstance(node.func, ast.Name):
                skip.add(id(node.func))
                if node.func.id == "Q" and len(node.args) == 1:
                    arg = node.args[0]
                    if isinstance(arg, ast.Constant) and isinstance(arg.value, str):
                        found.append(arg.value)
        elif isinstance(node, ast.Attribute) and isinstance(node.value, ast.Name):
            skip.add(id(node.value))
    for node in ast.walk(tree):
        if isinstance(node, ast.Name) and id(node) not in skip:
            found.append(node.id)
    return list(dict.fromkeys(found))


def _requests(mapping, known: set[str]) -> list:
    """Fields to resolve for the values of ``mapping``.

    Strings that are fields themselves (``"CD3D"``, ``"HLA-DRA"``, prefixed
    ``"gene:CD3D"``) are kept whole; any other string contributes the field
    identifiers it contains. Accessors (:class:`Ref`, :class:`ObsmRef`) are
    kept. Stages, callables and literals reference existing columns.
    """
    out: list = []
    for value in mapping.values():
        if isinstance(value, (Ref, ObsmRef)):
            out.append(value)
        elif isinstance(value, str):
            if isinstance(parse_token(value), (Ref, ObsmRef)) or value in known:
                out.append(value)
            else:
                out.extend(name for name in _identifiers(value) if name in known)
    return out


def _plain_mapping(adata, mapping):
    """Rewrite accessor / prefix-string aes values to their plain column names.

    plotnine evaluates aes values against the DataFrame, so ``"gene:CD3D"``,
    ``obsm("umap", 0)`` etc. must be swapped for the string they resolve to
    (``"CD3D"``, ``"UMAP_1"``) before the mapping reaches ``ggplot``. The
    DataFrame from :func:`resolve_frame` is keyed by exactly those names.
    """
    plain = aes()
    for key, value in mapping.items():
        if isinstance(value, (Ref, ObsmRef, str)):
            plain[key] = plain_name(adata, value)
        else:
            plain[key] = value
    return plain


def _facet_fields(facet) -> list[str]:
    names = list(getattr(facet, "vars", None) or [])
    names += list(getattr(facet, "rows", None) or []) + list(getattr(facet, "cols", None) or [])
    return [name for name in names if isinstance(name, str)]


# --------------------------------------------------------------------------- #
# The plot: a ggplot that remembers where its columns came from
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class _Source:
    """Where a :class:`GGAnnPlot` reads late fields from (immutable)."""

    adata: object
    layer: str | None
    use_raw: bool | None
    max_matrix_values: int | None
    spent: int  # matrix values materialized so far, against the budget
    obs_names: pd.Index  # row identity at creation, to detect an edited AnnData


class GGAnnPlot(ggplot):
    """A :class:`plotnine.ggplot` that resolves AnnData fields added with ``+``.

    Returned by :func:`gganndata`. ``plot + aes(color="louvain")``,
    ``plot + geom_point(aes(color="CD3D"))`` and ``plot + facet_wrap("~phase")``
    read the new field from the AnnData the plot was built from, within the
    original ``max_matrix_values`` budget. ``+`` still returns a new plot and the
    original is unchanged; the AnnData is never modified.

    Notes
    -----
    Late fields are read when they are added, not when the plot was created: an
    AnnData edited in between is not snapshotted (removing or reordering cells
    raises). A geom with its own ``data=`` is evaluated against that data and
    requests nothing.
    """

    _ggann_source: _Source | None = None

    def __deepcopy__(self, memo):
        result = super().__deepcopy__(memo)
        result.__dict__["_ggann_source"] = self.__dict__.get("_ggann_source")
        return result

    def __iadd__(self, other):
        source = self.__dict__.get("_ggann_source")
        if source is None:
            return super().__iadd__(other)
        return super().__iadd__(self._resolve_late(source, other))

    def _resolve_late(self, source: _Source, other):
        if isinstance(other, (list, tuple)):
            return [self._resolve_late(source, item) for item in other]
        adata = source.adata
        is_aes = isinstance(other, aes)
        is_layer = (
            isinstance(getattr(other, "mapping", None), aes) and getattr(other, "data", 0) is None
        )
        facet_fields = [] if is_aes or is_layer else _facet_fields(other)
        if not (is_aes or is_layer or facet_fields):
            return other  # themes, scales, labels, layers with their own data
        known = field_names(adata, layer=source.layer, use_raw=source.use_raw)
        if is_aes:
            requests = _requests(other, known)
        elif is_layer:
            requests = _requests(other.mapping, known)
        else:
            requests = [f for f in facet_fields if f in known]
        wanted = [r for r in requests if plain_name(adata, r) not in self.data.columns]
        if wanted:
            self._extend(source, wanted)
        if is_aes:
            return _plain_mapping(adata, other)
        if is_layer:
            layer = copy(other)  # the user's geom object is never edited
            layer.mapping = _plain_mapping(adata, other.mapping)
            return layer
        return other

    def _extend(self, source: _Source, requests: list) -> None:
        adata = source.adata
        if adata.n_obs != len(self.data) or not adata.obs_names.equals(source.obs_names):
            raise ValueError(
                "The AnnData changed since gganndata(): its cells are no longer the "
                "rows of this plot, so later fields cannot be added. Rebuild the plot."
            )
        remaining = (
            None if source.max_matrix_values is None else source.max_matrix_values - source.spent
        )
        cost: list[int] = []
        new = resolve_frame(
            adata,
            requests,
            layer=source.layer,
            use_raw=source.use_raw,
            max_matrix_values=remaining,
            _cost=cost,
        )
        new = new.loc[:, [c for c in new.columns if c not in self.data.columns]]
        data = self.data.copy(deep=False)
        for column in new.columns:  # by position: obs names may repeat
            data[column] = new[column].to_numpy()
        self.data = data
        self.__dict__["_ggann_source"] = replace(source, spent=source.spent + sum(cost))


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
        ``"UMAP_1"`` or ``"PC_5"`` (any column of any ``obsm`` array). Use
        :func:`ggann.gene` / :func:`ggann.obs` to force a source. A string may
        also be an expression over fields, e.g. ``"np.log1p(n_counts)"`` or
        ``"factor(louvain)"``; see Notes.
    layer
        Read expression from ``adata.layers[layer]`` instead of ``X``.
    use_raw
        Read expression from ``adata.raw``. Defaults to ``True`` when ``adata.raw``
        exists and no ``layer`` is given (scanpy's convention).
    max_matrix_values
        Maximum cumulative number of expression and ``obsm`` values that may be
        materialized for the whole plot, including fields added later with ``+``.
        The request is rejected before extraction when it exceeds this
        annplyr-compatible budget. Observation metadata does not count toward
        the budget. ``None`` disables the limit.
    add_theme
        Add :func:`theme_ggann` to the returned plot. With ``False`` the plot
        keeps plotnine's default ``theme_gray``; the return type is the same.

    Returns
    -------
    plotnine.ggplot
        A :class:`GGAnnPlot`. Compose it with any plotnine layer, e.g.
        ``gganndata(adata, aes("UMAP_1", "UMAP_2", color="louvain")) + geom_point()``.

    Raises
    ------
    KeyError
        If an explicitly selected observation, gene, embedding, or layer is absent.
    ValueError
        If ``layer`` and ``use_raw=True`` are combined, or the AnnData lost or
        reordered cells before a field was added with ``+``.
    annplyr.AnnplyrError
        If ``max_matrix_values`` is negative or the cumulative projection exceeds it.

    Notes
    -----
    Only referenced columns are projected through annplyr. The input ``adata``
    is not mutated.

    Expressions: field identifiers are found by parsing the string (never by
    evaluating it); plotnine then evaluates it against the plot data, in the
    namespace of the code that called :func:`gganndata`, plus ``factor`` and
    ``reorder``. A name that is neither a field nor defined there fails at draw
    time with plotnine's "could not evaluate" error. A gene whose name is not a
    Python identifier (``HLA-DRA``) is written ``gene("HLA-DRA")`` or, inside an
    expression, ``Q("HLA-DRA")``. Prefix tokens (``"gene:CD3D"``) work only as
    the whole value, not inside an expression.

    Examples
    --------
    >>> from plotnine import geom_point
    >>> p = gganndata(adata, aes("UMAP_1", "UMAP_2", color="cell_type"))
    >>> p = p + geom_point()
    """
    mapping = aes() if mapping is None else mapping
    cost: list[int] = []
    known = field_names(adata, layer=layer, use_raw=use_raw)
    df = resolve_frame(
        adata,
        _requests(mapping, known),
        layer=layer,
        use_raw=use_raw,
        max_matrix_values=max_matrix_values,
        _cost=cost,
    )
    plot = GGAnnPlot(df, _plain_mapping(adata, mapping))
    # evaluate expressions in the user's namespace, not in this module's
    caller = Environment.capture(1)
    plot.environment = Environment([*caller.namespaces, {"Q": Q}])  # Q("HLA-DRA") in expressions
    plot.__dict__["_ggann_source"] = _Source(
        adata=adata,
        layer=layer,
        use_raw=use_raw,
        max_matrix_values=max_matrix_values,
        spent=sum(cost),
        obs_names=adata.obs_names.copy(),
    )
    if add_theme:
        plot = plot + _family_theme("standard")
    return plot
