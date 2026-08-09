"""Gene-weighted kernel density on an embedding (pyNebulosa), rendered plotnine-native.

Standard feature plots are hurt by dropout: many cells read zero even for real
markers. pyNebulosa recovers the signal with a weighted 2D kernel density
estimate over the embedding. This module reuses pyNebulosa's ``calculate_density``
for the KDE but renders the result as an ordinary :class:`plotnine.ggplot`, so
density plots compose with ggann's theme / scales / facets -- unlike pyNebulosa's
own matplotlib ``plot_density``.

Expression weights are pulled through annplyr (via :func:`resolve_frame`), never
by indexing ``adata.X`` directly, keeping the one-extraction-layer contract.
"""

from __future__ import annotations

from typing import Sequence

import numpy as np
import pandas as pd
import plotnine_extra as pe
from plotnine import aes, geom_point, ggplot, labs, scale_color_cmap

from ._expression import ordered_unique
from ._resolve import embedding_coords, plain_name, resolve_frame
from .theme import theme_ggann

__all__ = ["plot_density"]


def _require_pynebulosa():
    try:
        from pynebulosa import calculate_density
    except ImportError as exc:  # pragma: no cover - exercised only without the dep
        raise ImportError(
            "plot_density requires pynebulosa; install with `pip install pynebulosa` "
            "(bundled in the ggann[density] extra)."
        ) from exc
    return calculate_density


def _minmax(values: np.ndarray) -> np.ndarray:
    """Scale to ``[0, 1]`` for display; monotonic, so per-cell ranking is preserved."""
    values = np.asarray(values, dtype=float)
    low, high = float(np.min(values)), float(np.max(values))
    if high <= low:
        return np.zeros_like(values)
    return (values - low) / (high - low)


def plot_density(
    adata,
    features: str | Sequence[str],
    *,
    joint: bool = False,
    basis: str = "umap",
    layer: str | None = None,
    use_raw: bool | None = None,
    method: str = "wkde",
    adjust: float = 1.0,
    size: float = 1.5,
    alpha: float = 0.9,
    cmap: str = "magma",
    ncol: int | None = None,
):
    """Gene-weighted density over an embedding, one faceted panel per feature.

    Reuses pyNebulosa's weighted KDE (``method="wkde"`` or ``"ks"``) to smooth
    sparse expression over the embedding, recovering marker signal lost to
    dropout. ``features`` may be genes or numeric obs columns.

    With ``joint=True`` and multiple features, a final panel shows the **joint
    density** -- the element-wise *product* of the individual densities (matching
    pyNebulosa / the R Nebulosa package), highlighting cells co-expressing all
    features.

    Each panel's colour is the density min-max scaled to ``[0, 1]`` so panels
    stay legible under plotnine's shared facet scale (raw density magnitudes
    differ by orders of magnitude across genes). The scaling is monotonic, so the
    within-panel ordering is unchanged; for raw-magnitude colourbars, compose
    separate ``plot_density`` calls with the re-exported ``Wrap`` / ``plot_layout``.

    Parameters
    ----------
    adata
        Annotated data matrix.
    features : str or sequence of str
        Genes or numeric observation columns.
    joint : bool
        Add the product-density panel for multiple features.
    basis : str
        Embedding basis.
    layer, use_raw : optional
        Mutually exclusive expression source.
    method : {"wkde", "ks"}
        pyNebulosa density estimator.
    adjust : float
        Bandwidth multiplier.
    size, alpha : float
        Point size and opacity.
    cmap : str
        Matplotlib colormap name.
    ncol : int, optional
        Facet columns.

    Returns
    -------
    plotnine.ggplot
        Composable faceted density plot.

    Raises
    ------
    ImportError
        If the ``density`` extra is unavailable.
    KeyError
        If the embedding or an explicit feature source is missing.
    TypeError
        If a resolved feature is non-numeric.
    ValueError
        If the method or expression-source selection is invalid.

    Notes
    -----
    Only selected features and two coordinates are projected. KDE cost scales with
    observations and features; the input is not mutated.

    Examples
    --------
    >>> p = plot_density(adata, ["CD3D", "NKG7"], joint=True)
    """
    calculate_density = _require_pynebulosa()

    if isinstance(features, str):
        features = [features]
    features = ordered_unique(features)
    if not features:
        raise ValueError("plot_density needs at least one feature.")

    coords = embedding_coords(adata, basis)
    if coords.shape[1] < 2:
        raise ValueError(f"Embedding '{basis}' has fewer than 2 dimensions.")
    xcol, ycol = coords.columns[:2]
    xy = coords.iloc[:, :2].to_numpy(dtype=float)

    values = resolve_frame(adata, features, layer=layer, use_raw=use_raw)
    names = [plain_name(adata, f) for f in features]
    missing = [f for f, n in zip(features, names) if n not in values.columns]
    if missing:
        raise KeyError(f"Could not resolve feature(s) {missing} as genes or obs columns.")
    non_numeric = [
        f for f, n in zip(features, names) if not pd.api.types.is_numeric_dtype(values[n])
    ]
    if non_numeric:
        raise TypeError(
            f"plot_density needs numeric features (genes or continuous obs); "
            f"{non_numeric} are non-numeric. Use plot_embedding for categorical colours."
        )

    panels = list(dict.fromkeys(names))
    raw_density = {
        name: np.asarray(
            calculate_density(
                np.asarray(values[name].to_numpy(), dtype=float),
                xy,
                method=method,
                adjust=adjust,
            ),
            dtype=float,
        )
        for name in panels
    }
    panel_density = dict(raw_density)
    if joint and len(panels) > 1:
        joint_label = " + ".join(panels)
        panel_density[joint_label] = np.prod([raw_density[n] for n in panels], axis=0)
        panels.append(joint_label)

    frames = [
        pd.DataFrame(
            {
                xcol: xy[:, 0],
                ycol: xy[:, 1],
                "density": _minmax(panel_density[label]),
                "feature": label,
            }
        )
        for label in panels
    ]
    long = pd.concat(frames, ignore_index=True)
    long["feature"] = pd.Categorical(long["feature"], categories=panels, ordered=True)
    long = long.sort_values(["feature", "density"])  # draw low density first, within each panel

    return (
        ggplot(long, aes(xcol, ycol, color="density"))
        + geom_point(size=size, alpha=alpha)
        + pe.facet_wrap("~feature", ncol=ncol)
        + scale_color_cmap(cmap_name=cmap)
        + labs(color="density\n(scaled)")
        + theme_ggann()
    )
