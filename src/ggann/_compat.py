"""Compatibility helpers for renamed keyword arguments."""

from __future__ import annotations

import warnings


def renamed_keyword(
    value,
    legacy_value,
    *,
    name: str,
    legacy_name: str,
    default,
):
    """Resolve a canonical keyword and one deprecated spelling."""
    if legacy_value is not None:
        if value is not None:
            raise TypeError(f"Pass only {name!r}; do not also pass {legacy_name!r}.")
        warnings.warn(
            f"{legacy_name!r} is deprecated; use {name!r} instead.",
            FutureWarning,
            stacklevel=3,
        )
        return legacy_value
    return default if value is None else value
