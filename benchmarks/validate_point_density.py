#!/usr/bin/env python3
"""Validate the binned point density against the exact KDE (decision 3, row P2).

Compares ``plot_embedding(method="binned")`` with the exact ``stat_pointdensity``
KDE (statsmodels, ``normal_reference`` bandwidth) at every cell, on scanpy's
bundled pbmc68k_reduced embedding and on a synthetic 20k-cell mixture.

Run: python benchmarks/validate_point_density.py
Prints one line per fixture: n, Spearman rank correlation, max normalized
error (max |binned - exact| / max exact), and both timings.
"""

from __future__ import annotations

import time

import numpy as np
import scanpy as sc
from plotnine.stats.density import kde
from scipy.stats import spearmanr

from ggann.plots import _binned_density, _kde_params


def exact_density(xy: np.ndarray) -> np.ndarray:
    return kde(xy, xy, "statsmodels-m", var_type="cc", **_kde_params()).ravel()


def synthetic(n: int = 20_000, seed: int = 0) -> np.ndarray:
    rng = np.random.default_rng(seed)
    centres = rng.uniform(-8, 8, size=(8, 2))
    which = rng.integers(0, len(centres), n)
    scale = rng.uniform(0.4, 1.2, size=len(centres))[which, None]
    return centres[which] + rng.normal(size=(n, 2)) * scale


fixtures = {
    "pbmc68k_reduced (UMAP)": np.asarray(sc.datasets.pbmc68k_reduced().obsm["X_umap"], float),
    "synthetic 20k mixture": synthetic(),
}
for name, xy in fixtures.items():
    t0 = time.perf_counter()
    exact = exact_density(xy)
    t_exact = time.perf_counter() - t0
    t0 = time.perf_counter()
    binned = _binned_density(xy)
    t_binned = time.perf_counter() - t0
    rho = spearmanr(exact, binned).statistic
    err = float(np.abs(binned - exact).max() / exact.max())
    print(
        f"{name}: n={len(xy)} spearman={rho:.4f} max_norm_err={err:.4f} "
        f"exact={t_exact:.2f}s binned={t_binned:.3f}s"
    )
