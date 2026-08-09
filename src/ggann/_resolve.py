"""Aesthetic resolution: turn names into a tidy per-cell DataFrame via annplyr.

Every plot in ggann is built from a plain :class:`pandas.DataFrame` that is
extracted from the ``AnnData`` *only* through the ``adata.ap`` (annplyr)
accessor -- no direct indexing into ``adata.X`` / ``adata.obs``. This module is
that single extraction layer; the plotting helpers sit on top of it.

Resolution order for a bare name (matching scanpy's precedence, least
surprising first):

1. a column of ``adata.obs``            (per-cell metadata)
2. a gene / feature in ``X`` or a layer (or ``adata.raw`` when ``use_raw``)
3. an embedding coordinate in ``obsm``  (e.g. ``"UMAP_1"``)

When a name is both an obs column and a gene, obs wins and a warning is
emitted. Use the :func:`gene` and :func:`obs` escapes to force a source.
"""

from __future__ import annotations

import operator
import re
import warnings
from dataclasses import dataclass
from typing import Iterable

import pandas as pd

from ._expression import (
    _expression_kwargs,
    expression_frame,
    resolve_source,
    source_var_names,
    validate_materialization_budget,
)
from ._expression import (
    densify_frame as _densify,
)
from ._expression import (
    source_label as _source_label,
)

__all__ = ["gene", "obs", "obsm", "embedding_coords"]

_OBSM_TOKEN = re.compile(r"^(?P<basis>.+?)\[(?P<idx>\d+)\]$")


@dataclass(frozen=True)
class Ref:
    """An explicit reference to a data source, produced by :func:`gene`/:func:`obs`."""

    name: str
    source: str  # "gene" | "obs"
    layer: str | None = None
    use_raw: bool | None = None

    def __str__(self) -> str:  # so it can be dropped straight into aes()
        return self.name


def gene(name: str, *, layer: str | None = None, use_raw: bool | None = None) -> Ref:
    """Force ``name`` to resolve as a gene/feature (expression), never obs.

    Optionally pin *this gene's* expression matrix, independent of the plot-wide
    ``layer`` / ``use_raw``. This is how you mix sources in one plot::

        aes(color=gene("CD3D", layer="logcounts"),
            size=gene("CD8A", layer="counts"))

    With neither ``layer`` nor ``use_raw`` given, the gene inherits the plot-wide
    source (which itself defaults to ``adata.raw`` when present).

    Parameters
    ----------
    name : str
        Variable name in the selected expression source.
    layer : str, optional
        Named ``adata.layers`` source for this reference.
    use_raw : bool, optional
        Select ``adata.raw`` (true) or ``adata.X`` (false) for this reference.

    Returns
    -------
    Ref
        An immutable source-qualified aesthetic reference.

    Raises
    ------
    ValueError
        When resolution later combines ``layer`` with ``use_raw=True``.

    Notes
    -----
    Construction does not inspect or mutate an AnnData object.

    Examples
    --------
    >>> gene("CD3D", layer="log1p")
    """
    return Ref(str(name), "gene", layer=layer, use_raw=use_raw)


def obs(name: str) -> Ref:
    """Force ``name`` to resolve as an ``adata.obs`` column, never a gene.

    Parameters
    ----------
    name : str
        Observation-column name.

    Returns
    -------
    Ref
        An immutable observation reference.

    Raises
    ------
    KeyError
        When the reference is later resolved against data without ``name``.

    Notes
    -----
    Construction does not inspect or mutate an AnnData object.

    Examples
    --------
    >>> obs("cell_type")
    """
    return Ref(str(name), "obs")


@dataclass(frozen=True)
class ObsmRef:
    """A reference to one embedding coordinate, produced by :func:`obsm`."""

    basis: str
    index: int

    def __str__(self) -> str:
        return f"{self.basis}[{self.index}]"


def obsm(basis: str, index: int) -> ObsmRef:
    """Force a specific zero-based embedding coordinate.

    Parameters
    ----------
    basis : str
        User-facing basis or exact ``adata.obsm`` key.
    index : int
        Zero-based coordinate index.

    Returns
    -------
    ObsmRef
        An immutable embedding-coordinate reference.

    Raises
    ------
    TypeError
        If ``index`` is not an integer.
    KeyError
        When the basis is later resolved but is absent.
    IndexError
        When the coordinate is later resolved but is out of range.

    Notes
    -----
    Construction does not read or mutate an AnnData object.

    Examples
    --------
    >>> obsm("umap", 0)
    """
    if isinstance(index, bool):
        raise TypeError(f"obsm index must be an integer, got {index!r}.")
    try:
        index = operator.index(index)
    except TypeError as error:
        raise TypeError(f"obsm index must be an integer, got {index!r}.") from error
    return ObsmRef(str(basis), index)


# --------------------------------------------------------------------------- #
# Prefix-string parsing: "obs:phase", "gene:CD3D@logcounts", "obsm:umap[0]"
# --------------------------------------------------------------------------- #
def parse_token(item):
    """Normalise one aes value into a :class:`Ref` / :class:`ObsmRef` / bare string.

    Bare strings (no recognised prefix) are returned unchanged for
    auto-resolution. Recognised strict prefixes:

    * ``obs:<col>``               -- an ``adata.obs`` column
    * ``gene:<name>``             -- expression, using the plot-wide source
    * ``gene:<name>@<layer>``     -- expression from ``adata.layers[<layer>]``
    * ``gene:<name>@raw``/``@X``  -- expression from ``adata.raw`` / ``adata.X``
    * ``obsm:<basis>[<i>]``       -- coordinate ``i`` (0-based) of an embedding

    Non-string, non-Ref values are returned unchanged.
    """
    if isinstance(item, (Ref, ObsmRef)):
        return item
    if not isinstance(item, str) or ":" not in item:
        return item
    prefix, _, body = item.partition(":")
    kind = prefix.strip().lower()
    if kind == "obs":
        return Ref(body.strip(), "obs")
    if kind == "gene":
        name, at, source = body.partition("@")
        name, source = name.strip(), source.strip()
        if not at:
            return Ref(name, "gene")
        low = source.lower()
        if low == "raw":
            return Ref(name, "gene", use_raw=True)
        if low in ("x", ".x"):
            return Ref(name, "gene", use_raw=False)
        return Ref(name, "gene", layer=source)
    if kind == "obsm":
        m = _OBSM_TOKEN.match(body.strip())
        if not m:
            raise ValueError(f"obsm token {item!r} must look like 'obsm:umap[0]' (0-based index).")
        return ObsmRef(m.group("basis"), int(m.group("idx")))
    return item  # unrecognised prefix -> leave for plotnine / auto-resolution


def plain_name(adata, item):
    """The plain DataFrame column name a token resolves to (for aes rewriting)."""
    tok = parse_token(item)
    if isinstance(tok, ObsmRef):
        key = embedding_key(adata, tok.basis)
        return f"{_embedding_prefix(key)}_{tok.index + 1}"
    if isinstance(tok, Ref):
        return tok.name
    return tok


def _other_gene_universe(adata, kind: str) -> pd.Index:
    if kind == "raw":
        return adata.var_names
    return adata.raw.var_names if adata.raw is not None else pd.Index([])


# --------------------------------------------------------------------------- #
# Embeddings
# --------------------------------------------------------------------------- #
def _embedding_prefix(key: str) -> str:
    """``X_umap`` -> ``UMAP``, ``X_pca`` -> ``PC`` (Seurat-style coordinate names)."""
    base = key[2:] if key.startswith("X_") else key
    base = base.upper()
    return {"PCA": "PC"}.get(base, base)


def embedding_key(adata, basis: str) -> str:
    """Resolve a user-facing basis name (``"umap"``, ``"X_umap"``, ``"UMAP"``) to an obsm key."""
    candidates = [basis, f"X_{basis}", basis.lower(), f"X_{basis.lower()}"]
    for cand in candidates:
        if cand in adata.obsm:
            return cand
    lower = {k.lower(): k for k in adata.obsm}
    for cand in (basis.lower(), f"x_{basis.lower()}"):
        if cand in lower:
            return lower[cand]
    raise KeyError(f"No embedding '{basis}' in adata.obsm (available: {list(adata.obsm)})")


def embedding_coords(
    adata,
    basis: str,
    n: int = 2,
    *,
    max_matrix_values: int | None = None,
) -> pd.DataFrame:
    """Return the first ``n`` coordinates of an embedding as a tidy DataFrame.

    Columns are named ``<PREFIX>_<i>`` (e.g. ``UMAP_1``, ``UMAP_2``) to match the
    conventions used by ``plotnine_extra.DimPlot``.

    Parameters
    ----------
    adata
        Annotated data matrix containing the embedding.
    basis : str
        User-facing basis or exact ``adata.obsm`` key.
    n : int, default=2
        Maximum number of leading coordinates to project.
    max_matrix_values : int, optional
        Maximum number of embedding values allowed before extraction.

    Returns
    -------
    pandas.DataFrame
        Projected coordinates in observation order, including duplicate names.

    Raises
    ------
    KeyError
        If ``basis`` is absent.
    ValueError
        If ``n`` is negative.
    annplyr.AnnplyrError
        If the materialization budget is invalid or exceeded.

    Notes
    -----
    annplyr projects only the requested coordinates. ``adata`` is not mutated.

    Examples
    --------
    >>> coords = embedding_coords(adata, "umap", n=2)
    """
    if n < 0:
        raise ValueError(f"n must be non-negative, got {n}.")
    key = embedding_key(adata, basis)
    width = adata.obsm[key].shape[1]
    indices = range(min(n, width))
    validate_materialization_budget(
        max_matrix_values,
        adata.n_obs * min(n, width),
        context="embedding_coords",
    )
    return _embedding_frame(
        adata,
        key,
        indices,
        max_matrix_values=max_matrix_values,
    )


def _embedding_frame(
    adata,
    key: str,
    indices: Iterable[int],
    *,
    max_matrix_values: int | None = None,
) -> pd.DataFrame:
    """Extract selected embedding coordinates through annplyr."""
    width = adata.obsm[key].shape[1]
    indices = list(dict.fromkeys(indices))
    invalid = [index for index in indices if index < 0 or index >= width]
    if invalid:
        raise IndexError(
            f"Embedding {key!r} has {width} coordinates; requested index {invalid[0]}."
        )
    if not indices:
        return pd.DataFrame(index=adata.obs_names)
    selectors = _embedding_selectors(adata, key, indices)
    frame = adata.ap.to_df(
        obsm={key: selectors},
        max_matrix_values=max_matrix_values,
    )
    prefix = _embedding_prefix(key)
    frame.columns = [f"{prefix}_{index + 1}" for index in indices]
    return _densify(frame)


def _embedding_selectors(adata, key: str, indices: Iterable[int]) -> list[str]:
    """Return annplyr selectors for positional embedding coordinates."""
    indices = list(indices)
    value = adata.obsm[key]
    return (
        [str(value.columns[index]) for index in indices]
        if isinstance(value, pd.DataFrame)
        else [str(index) for index in indices]
    )


def _all_embedding_coords(adata, n: int = 2) -> dict[str, tuple[str, int]]:
    """Map every embedding coordinate name -> (obsm key, column index)."""
    out: dict[str, tuple[str, int]] = {}
    for key, arr in adata.obsm.items():
        if getattr(arr, "ndim", 0) != 2:
            continue
        prefix = _embedding_prefix(key)
        for i in range(min(n, arr.shape[1])):
            out.setdefault(f"{prefix}_{i + 1}", (key, i))
    return out


# --------------------------------------------------------------------------- #
# Frame resolution
# --------------------------------------------------------------------------- #
def resolve_frame(
    adata,
    names: Iterable,
    *,
    layer: str | None = None,
    use_raw: bool | None = None,
    warn_collisions: bool = True,
    max_matrix_values: int | None = None,
) -> pd.DataFrame:
    """Build a per-cell DataFrame with one column per requested ``name``.

    ``names`` may contain plain strings or :class:`Ref` objects (from
    :func:`gene`/:func:`obs`). Names that resolve to nothing are silently
    skipped -- they are assumed to be computed aesthetics that plotnine will
    handle from already-present columns. A name that is a gene only in the
    *inactive* matrix (e.g. in ``.X`` while reading from ``.raw``) is skipped
    *with a warning*, so the mistake is visible.
    """
    default_kind, default_layer = resolve_source(adata, layer, use_raw)
    default_set = source_var_names(adata, default_kind)
    other_set = _other_gene_universe(adata, default_kind)
    obs_cols = set(adata.obs.columns)
    emb = _all_embedding_coords(adata)

    obs_names: list[str] = []
    gene_specs: list[tuple[str, str, str | None]] = []  # (name, kind, layer)
    obsm_specs: list[tuple[str, str, int]] = []  # (col_name, obsm_key, index)
    order: list[str] = []

    def _ref_source(ref: Ref) -> tuple[str, str | None]:
        # an explicit per-gene layer/use_raw fully determines the source;
        # otherwise inherit the plot-wide default.
        if ref.layer is not None:
            return resolve_source(adata, ref.layer, ref.use_raw)
        if ref.use_raw is not None:
            return resolve_source(adata, None, ref.use_raw)
        return default_kind, default_layer

    for raw_item in names:
        tok = parse_token(raw_item)

        if isinstance(tok, ObsmRef):
            key = embedding_key(adata, tok.basis)
            name = f"{_embedding_prefix(key)}_{tok.index + 1}"
            if name in order:
                continue
            obsm_specs.append((name, key, tok.index))
            order.append(name)
            continue

        if isinstance(tok, Ref):
            name = tok.name
            if name in order:
                continue
            if tok.source == "obs":
                _dispatch(name, "obs", obs_cols, obs_names)
            else:  # forced gene
                k, lyr = _ref_source(tok)
                if name not in source_var_names(adata, k):
                    raise KeyError(f"gene('{name}') not found in {_source_label(k, lyr)}.")
                gene_specs.append((name, k, lyr))
            order.append(name)
            continue

        # bare string -> auto-resolve (obs -> gene -> obsm coordinate)
        name = str(tok)
        if name in order:
            continue
        if name in obs_cols:
            if warn_collisions and name in default_set:
                warnings.warn(
                    f"'{name}' is both an obs column and a gene; using obs. "
                    f"Use 'gene:{name}' to plot expression.",
                    stacklevel=2,
                )
            obs_names.append(name)
        elif name in default_set:
            gene_specs.append((name, default_kind, default_layer))
        elif name in emb:
            key, idx = emb[name]
            obsm_specs.append((name, key, idx))
        elif name in other_set:
            warnings.warn(
                f"'{name}' is a gene in "
                f"{_source_label('raw' if default_kind != 'raw' else 'x', default_layer)} "
                f"but not in the active expression source "
                f"({_source_label(default_kind, default_layer)}); "
                f"pass use_raw={default_kind != 'raw'} to plot it.",
                stacklevel=2,
            )
            continue
        else:
            continue  # computed / literal aesthetic -- leave to plotnine
        order.append(name)

    # group genes by their (matrix, layer) so mixed per-aesthetic sources each
    # get one extraction pass
    by_source: dict[tuple[str, str | None], list[str]] = {}
    for name, k, lyr in gene_specs:
        by_source.setdefault((k, lyr), []).append(name)

    # A single aesthetic may span X, raw, several layers, and several obsm keys.
    # annplyr preflights cumulative budgets within one call; ggann preflights the
    # complete multi-call request here before the first extraction.
    projected_columns = sum(len(gnames) for gnames in by_source.values())
    by_key: dict[str, list[tuple[str, int]]] = {}
    for name, key, idx in obsm_specs:
        by_key.setdefault(key, []).append((name, idx))
    projected_columns += sum(len(cols) for cols in by_key.values())
    validate_materialization_budget(
        max_matrix_values,
        adata.n_obs * projected_columns,
        context="resolve_frame",
    )

    frames = []

    # annplyr can preflight and execute obs, one X/layer projection, raw, and
    # multiple obsm projections together.  Coalesce that compatible subset so
    # common grammar/embedding plots cross the public extraction boundary once.
    source_items = list(by_source.items())
    fused_source_items: list[tuple[tuple[str, str | None], list[str]]] = []
    raw_item = next((item for item in source_items if item[0][0] == "raw"), None)
    matrix_item = next((item for item in source_items if item[0][0] != "raw"), None)
    if matrix_item is not None:
        fused_source_items.append(matrix_item)
    if raw_item is not None:
        fused_source_items.append(raw_item)

    obsm_requests: dict[str, list[str]] = {}
    for key, cols in by_key.items():
        indices = [index for _, index in cols]
        width = adata.obsm[key].shape[1]
        invalid = [index for index in indices if index < 0 or index >= width]
        if invalid:
            raise IndexError(
                f"Embedding {key!r} has {width} coordinates; requested index {invalid[0]}."
            )
        obsm_requests[key] = _embedding_selectors(adata, key, indices)

    exported_names = list(obs_names)
    for (k, _lyr), gnames in fused_source_items:
        exported_names.extend(f"raw_{name}" if k == "raw" else name for name in gnames)
    for key, selectors in obsm_requests.items():
        exported_names.extend(f"{key}_{selector}" for selector in selectors)
    can_fuse = bool(obs_names or fused_source_items or obsm_requests) and len(
        exported_names
    ) == len(set(exported_names))

    fused_keys: set[tuple[str, str | None]] = set()
    if can_fuse:
        call_kwargs = {"obs": obs_names or None, "obsm": obsm_requests or None}
        layout: list[tuple[list[str], int]] = []
        if obs_names:
            layout.append((obs_names, len(obs_names)))
        # annplyr emits X/layer before raw regardless of keyword insertion order.
        for expected_kind in ("x_or_layer", "raw"):
            for (k, lyr), gnames in fused_source_items:
                if (expected_kind == "raw") != (k == "raw"):
                    continue
                call_kwargs.update(_expression_kwargs(k, lyr, gnames))
                fused_keys.add((k, lyr))
                layout.append((gnames, len(gnames)))
        for key, cols in by_key.items():
            layout.append(([name for name, _ in cols], len(cols)))

        combined = _densify(
            adata.ap.to_df(
                **call_kwargs,
                max_matrix_values=max_matrix_values,
            )
        )
        cursor = 0
        for names, width in layout:
            piece = combined.iloc[:, cursor : cursor + width].copy(deep=False)
            piece.columns = names
            frames.append(piece)
            cursor += width
    else:
        if obs_names:
            frames.append(_densify(adata.ap.to_df(obs=obs_names)))
        for key, cols in by_key.items():
            indices = [index for _, index in cols]
            coords = _embedding_frame(
                adata,
                key,
                indices,
                max_matrix_values=max_matrix_values,
            )
            coords.columns = [name for name, _ in cols]
            frames.append(coords)

    for (k, lyr), gnames in source_items:
        if (k, lyr) in fused_keys:
            continue
        frame, gnames = expression_frame(
            adata,
            gnames,
            kind=k,
            layer=lyr,
            max_matrix_values=max_matrix_values,
        )
        frame = _densify(frame)
        frames.append(frame)

    if not frames:
        return pd.DataFrame(index=adata.obs_names)

    # Assemble by position.  Label-aligned concat is ambiguous for duplicate
    # observation names even though each annplyr result preserves row order.
    out = pd.concat([frame.reset_index(drop=True) for frame in frames], axis=1)
    out.index = adata.obs_names.copy()
    ordered = [c for c in order if c in out.columns]
    return out[ordered + [c for c in out.columns if c not in ordered]]


def _dispatch(name, kind, universe, bucket):
    if name not in universe:
        raise KeyError(f"{kind}('{name}') not found in adata.")
    bucket.append(name)
