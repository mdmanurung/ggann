"""Keyword helpers: accept Scanpy spellings of ggann arguments."""

from __future__ import annotations

import functools
import inspect
import os

# Attribute warnings to the first frame outside ggann, at any call depth.
SKIP_GGANN_FRAMES = (os.path.dirname(__file__),)


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


# Scanpy spelling -> candidate ggann parameter names, first one the helper has.
_SCANPY_ALIASES = {
    "groupby": ("group_by",),
    "var_names": ("genes", "features"),
    "groups": ("group",),
    "cmap": ("color_map",),
    "color_map": ("cmap",),
}


def accepts_scanpy_names(func):
    """Accept Scanpy's spelling of a keyword as a silent alias.

    ``groupby`` maps to ``group_by``, ``var_names`` to ``genes`` or ``features``,
    and ``cmap`` / ``color_map`` to whichever of the two the helper spells.
    An alias is added only when the helper lacks that name and has the target.
    Both spellings with different values raise ``TypeError``.
    """
    params = inspect.signature(func).parameters
    aliases = {
        old: next(t for t in targets if t in params)
        for old, targets in _SCANPY_ALIASES.items()
        if old not in params and any(t in params for t in targets)
    }

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        for old, new in aliases.items():
            if old in kwargs:
                kwargs[new] = _resolve_alias(kwargs.get(new), kwargs.pop(old), new, old)
        return func(*args, **kwargs)

    return wrapper
