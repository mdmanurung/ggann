"""Helpers the gg* constructors in ``ggann.shapes`` share."""

from __future__ import annotations

import warnings

import pandas as pd
from plotnine import scale_color_manual, scale_fill_manual

from .._compat import SKIP_GGANN_FRAMES
from .._grouping import _order_groups, _resolve_group_order


def prepare_levels(
    df: pd.DataFrame, column: str, order, *, keep: list[str]
) -> tuple[pd.DataFrame, list]:
    """Copy ``keep`` columns of ``df`` and make ``column`` an ordered categorical.

    Order is the declared categorical order, else sorted, as in the rest of
    ggann; ``order`` overrides it. Returns the frame and the levels in use.
    """
    missing = [c for c in keep if c not in df.columns]
    if missing:
        raise KeyError(f"Column(s) not found in the data frame: {missing}.")
    out = df[list(dict.fromkeys(keep))].copy()
    out = _order_groups(out, column, _resolve_group_order(out, column, order))
    return out, list(out[column].cat.categories)


def drop_missing_rows(df: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    """Drop rows with a missing value in ``columns``, warning with the count."""
    columns = list(dict.fromkeys(c for c in columns if c is not None))
    missing = df[columns].isna().any(axis=1).to_numpy()
    n = int(missing.sum())
    if not n:
        return df
    names = ", ".join(repr(c) for c in columns)
    warnings.warn(
        f"dropped {n} of {len(df)} rows with missing {names}",
        UserWarning,
        skip_file_prefixes=SKIP_GGANN_FRAMES,
    )
    return df.loc[~missing]


def palette_colours(palette, levels: list) -> dict:
    """Map each level to a colour from a dict, a list or a colormap name."""
    from ..plots import _palette_colours

    if isinstance(palette, dict):
        lacking = [c for c in levels if c not in palette]
        if lacking:
            raise ValueError(f"palette has no colour for {lacking}.")
        return {c: palette[c] for c in levels}
    colours = _palette_colours(palette, len(levels))
    if len(colours) < len(levels):
        raise ValueError(f"palette has {len(colours)} colours for {len(levels)} levels.")
    return dict(zip(levels, colours))


def palette_scales(palette, levels: list, aesthetics: tuple[str, ...]) -> list:
    """Manual colour or fill scales for ``palette``; empty when it is ``None``."""
    if palette is None:
        return []
    values = palette_colours(palette, levels)
    makers = {"color": scale_color_manual, "fill": scale_fill_manual}
    return [makers[a](values=values) for a in aesthetics]
