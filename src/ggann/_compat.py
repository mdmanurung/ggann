"""Compatibility helpers for renamed keyword arguments."""

from __future__ import annotations

import functools
import os
import warnings

# Attribute warnings to the first frame outside ggann, at any call depth.
SKIP_GGANN_FRAMES = (os.path.dirname(__file__),)


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
            skip_file_prefixes=SKIP_GGANN_FRAMES,
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


def accepts_groupby(func):
    """Accept Scanpy's ``groupby=`` as a silent alias of ``group_by=``.

    Both spellings with different values raise ``TypeError``.
    """

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        if "groupby" in kwargs:
            kwargs["group_by"] = _resolve_alias(
                kwargs.get("group_by"), kwargs.pop("groupby"), "group_by", "groupby"
            )
        return func(*args, **kwargs)

    return wrapper
