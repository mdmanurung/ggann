"""Categorical ordering and deterministic observation sampling."""

from __future__ import annotations

import warnings
from collections.abc import Iterable

import numpy as np
import pandas as pd

from ._compat import SKIP_GGANN_FRAMES


def _group_categories(adata, group_by: str) -> list | None:
    """Return the declared order of a categorical observation column."""
    col = adata.obs[group_by]
    if isinstance(col.dtype, pd.CategoricalDtype):
        return list(col.cat.categories)
    return None


def _drop_missing_groups(adata, *fields):
    """Drop cells whose grouping field is missing, warning with the count.

    ``fields`` that are ``None`` or not ``obs`` columns are ignored, so callers
    can pass ``group_by`` and ``split_by`` unconditionally.
    """
    fields = [f for f in dict.fromkeys(fields) if isinstance(f, str) and f in adata.obs]
    if not fields:
        return adata
    missing = adata.obs[fields].isna().any(axis=1).to_numpy()
    n = int(missing.sum())
    if not n:
        return adata
    names = ", ".join(repr(f) for f in fields)
    warnings.warn(
        f"dropped {n} of {adata.n_obs} cells with missing {names}",
        UserWarning,
        skip_file_prefixes=SKIP_GGANN_FRAMES,
    )
    return adata[~missing]


def _check_group_cardinality(adata, field: str | None, limit: int, *, action: str) -> int:
    """Count the non-missing levels of ``field``; ``action`` is 'raise' or 'warn'.

    Guards helpers whose cost or readability grows with the number of groups
    (a numeric column such as a cell index has one level per cell).
    """
    if field is None or field not in adata.obs:
        return 0
    col = adata.obs[field]
    n = len(col.cat.categories) if isinstance(col.dtype, pd.CategoricalDtype) else col.nunique()
    n = int(min(n, col.nunique()))
    if n <= limit:
        return n
    kind = "numeric " if pd.api.types.is_numeric_dtype(col) else ""
    message = (
        f"{kind}group field {field!r} has {n} levels (limit {limit}); "
        "bin or convert it to a categorical with fewer levels first"
    )
    if action == "raise":
        raise ValueError(message)
    warnings.warn(message, UserWarning, skip_file_prefixes=SKIP_GGANN_FRAMES)
    return n


def _resolve_group_order(
    frame: pd.DataFrame,
    group_by: str,
    categories_order: Iterable | None,
) -> list:
    """Resolve and validate the categories used for one discrete variable."""
    if categories_order is not None:
        categories = list(dict.fromkeys(categories_order))
    else:
        col = pd.Series(frame[group_by], index=frame.index)
        if isinstance(col.dtype, pd.CategoricalDtype):
            categories = list(col.cat.categories)
        else:
            observed = list(pd.unique(col.dropna().to_numpy()))
            try:
                categories = sorted(observed)
            except TypeError:
                categories = observed

    present = list(pd.unique(pd.Series(frame[group_by]).dropna().astype(object)))
    category_set = set(categories)
    missing = [value for value in present if value not in category_set]
    if missing:
        raise ValueError(f"categories_order is missing groups present in the data: {missing}")
    return categories


def _order_groups(
    frame: pd.DataFrame,
    group_by: str,
    categories_order: Iterable | None,
) -> pd.DataFrame:
    """Apply a validated, unordered categorical axis to ``frame`` in place."""
    categories = _resolve_group_order(frame, group_by, categories_order)
    frame[group_by] = pd.Categorical(frame[group_by], categories=categories, ordered=False)
    frame[group_by] = frame[group_by].cat.remove_unused_categories()
    return frame


def _downsample_cells(
    adata,
    group_by: str | None,
    n: int | None,
    *,
    random_state: int | None = 0,
):
    """Return an observation view capped at ``n`` cells total or per group."""
    if n is None:
        return adata
    if isinstance(n, (bool, np.bool_)) or not isinstance(n, (int, np.integer)) or n < 1:
        raise ValueError(f"downsample must be a positive integer, got {n}.")
    if adata.n_obs <= n:
        return adata

    rng = np.random.RandomState(random_state)
    if group_by is None:
        keep = rng.choice(adata.n_obs, int(n), replace=False)
        return adata[np.sort(keep)]

    col = adata.obs[group_by]
    if isinstance(col.dtype, pd.CategoricalDtype):
        groups = list(col.cat.categories)
    else:
        groups = list(pd.unique(col.dropna()))

    parts: list[np.ndarray] = []
    sampled = False
    for group in groups:
        idx = np.flatnonzero(col.eq(group).fillna(False).to_numpy())
        if len(idx) > n:
            idx = rng.choice(idx, int(n), replace=False)
            sampled = True
        parts.append(idx)

    missing_idx = np.flatnonzero(col.isna().to_numpy())
    if len(missing_idx):
        if len(missing_idx) > n:
            missing_idx = rng.choice(missing_idx, int(n), replace=False)
            sampled = True
        parts.append(missing_idx)

    if not sampled:
        return adata
    keep = np.concatenate(parts) if parts else np.arange(adata.n_obs)
    return adata[np.sort(keep)]
