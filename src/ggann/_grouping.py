"""Categorical ordering and deterministic observation sampling."""

from __future__ import annotations

from collections.abc import Iterable

import numpy as np
import pandas as pd


def _group_categories(adata, group_by: str) -> list | None:
    """Return the declared order of a categorical observation column."""
    col = adata.obs[group_by]
    if isinstance(col.dtype, pd.CategoricalDtype):
        return list(col.cat.categories)
    return None


def _resolve_group_order(
    frame: pd.DataFrame,
    group_by: str,
    categories_order: Iterable | None,
) -> list:
    """Resolve and validate the categories used for one discrete variable."""
    if categories_order is not None:
        categories = list(dict.fromkeys(categories_order))
    else:
        col = frame[group_by]
        if isinstance(col.dtype, pd.CategoricalDtype):
            categories = list(col.cat.categories)
        else:
            observed = list(pd.unique(col.dropna()))
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
