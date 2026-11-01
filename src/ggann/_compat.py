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


def _resolve_alias(new_value, old_value, new_name: str, old_name: str):
    """Return the one value given for an alias pair; error if both differ."""
    if new_value is None:
        return old_value
    if old_value is not None and old_value != new_value:
        raise TypeError(
            f"{new_name!r} and {old_name!r} are aliases but were given different values "
            f"({new_value!r} vs {old_value!r}); pass only {new_name!r}."
        )
    return new_value
