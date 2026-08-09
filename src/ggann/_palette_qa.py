"""Deterministic colour-accessibility checks used by publication QA tests."""

from __future__ import annotations

import itertools
import math
from collections.abc import Sequence

import numpy as np
from matplotlib.colors import to_rgb

_CVD_MATRICES = {
    "protanopia": np.array(
        [
            [0.152286, 1.052583, -0.204868],
            [0.114503, 0.786281, 0.099216],
            [-0.003882, -0.048116, 1.051998],
        ]
    ),
    "deuteranopia": np.array(
        [
            [0.367322, 0.860646, -0.227968],
            [0.280085, 0.672501, 0.047413],
            [-0.011820, 0.042940, 0.968881],
        ]
    ),
    "tritanopia": np.array(
        [
            [1.255528, -0.076749, -0.178779],
            [-0.078411, 0.930809, 0.147602],
            [0.004733, 0.691367, 0.303900],
        ]
    ),
}


def _srgb_to_linear(rgb: np.ndarray) -> np.ndarray:
    return np.where(rgb <= 0.04045, rgb / 12.92, ((rgb + 0.055) / 1.055) ** 2.4)


def _linear_to_srgb(rgb: np.ndarray) -> np.ndarray:
    return np.where(
        rgb <= 0.0031308,
        12.92 * rgb,
        1.055 * np.power(np.clip(rgb, 0, None), 1 / 2.4) - 0.055,
    )


def simulate_cvd(rgb: np.ndarray, kind: str) -> np.ndarray:
    """Apply the Machado severity-100 linear-RGB approximation."""
    if kind not in _CVD_MATRICES:
        raise ValueError(f"Unknown colour-vision simulation {kind!r}.")
    linear = _srgb_to_linear(np.asarray(rgb, dtype=float))
    simulated = linear @ _CVD_MATRICES[kind].T
    return np.clip(_linear_to_srgb(np.clip(simulated, 0, 1)), 0, 1)


def _rgb_to_lab(rgb: np.ndarray) -> np.ndarray:
    linear = _srgb_to_linear(np.asarray(rgb, dtype=float))
    xyz = linear @ np.array(
        [
            [0.4124564, 0.2126729, 0.0193339],
            [0.3575761, 0.7151522, 0.1191920],
            [0.1804375, 0.0721750, 0.9503041],
        ]
    )
    xyz /= np.array([0.95047, 1.0, 1.08883])
    delta = 6 / 29
    f = np.where(xyz > delta**3, np.cbrt(xyz), xyz / (3 * delta**2) + 4 / 29)
    return np.column_stack(
        [116 * f[:, 1] - 16, 500 * (f[:, 0] - f[:, 1]), 200 * (f[:, 1] - f[:, 2])]
    )


def ciede2000(first: np.ndarray, second: np.ndarray) -> float:
    """CIEDE2000 colour difference for two CIELAB triples."""
    l1, a1, b1 = (float(value) for value in first)
    l2, a2, b2 = (float(value) for value in second)
    c1 = math.hypot(a1, b1)
    c2 = math.hypot(a2, b2)
    c_bar = (c1 + c2) / 2
    g = 0.5 * (1 - math.sqrt(c_bar**7 / (c_bar**7 + 25**7)))
    a1p, a2p = (1 + g) * a1, (1 + g) * a2
    c1p, c2p = math.hypot(a1p, b1), math.hypot(a2p, b2)

    def hue(a: float, b: float) -> float:
        angle = math.degrees(math.atan2(b, a))
        return angle + 360 if angle < 0 else angle

    h1p, h2p = hue(a1p, b1), hue(a2p, b2)
    dl = l2 - l1
    dc = c2p - c1p
    dh_angle = h2p - h1p
    if c1p * c2p == 0:
        dh_angle = 0
    elif dh_angle > 180:
        dh_angle -= 360
    elif dh_angle < -180:
        dh_angle += 360
    dh = 2 * math.sqrt(c1p * c2p) * math.sin(math.radians(dh_angle / 2))
    l_bar = (l1 + l2) / 2
    cp_bar = (c1p + c2p) / 2
    if c1p * c2p == 0:
        hp_bar = h1p + h2p
    elif abs(h1p - h2p) <= 180:
        hp_bar = (h1p + h2p) / 2
    elif h1p + h2p < 360:
        hp_bar = (h1p + h2p + 360) / 2
    else:
        hp_bar = (h1p + h2p - 360) / 2
    t = (
        1
        - 0.17 * math.cos(math.radians(hp_bar - 30))
        + 0.24 * math.cos(math.radians(2 * hp_bar))
        + 0.32 * math.cos(math.radians(3 * hp_bar + 6))
        - 0.20 * math.cos(math.radians(4 * hp_bar - 63))
    )
    sl = 1 + 0.015 * (l_bar - 50) ** 2 / math.sqrt(20 + (l_bar - 50) ** 2)
    sc = 1 + 0.045 * cp_bar
    sh = 1 + 0.015 * cp_bar * t
    delta_theta = 30 * math.exp(-(((hp_bar - 275) / 25) ** 2))
    rc = 2 * math.sqrt(cp_bar**7 / (cp_bar**7 + 25**7))
    rt = -rc * math.sin(math.radians(2 * delta_theta))
    return math.sqrt((dl / sl) ** 2 + (dc / sc) ** 2 + (dh / sh) ** 2 + rt * (dc / sc) * (dh / sh))


def palette_accessibility_report(colours: Sequence[str]) -> dict[str, object]:
    """Return grayscale luminance and minimum pairwise CIEDE2000 separations."""
    rgb = np.array([to_rgb(colour) for colour in colours], dtype=float)
    simulations = {"normal": rgb}
    simulations.update({name: simulate_cvd(rgb, name) for name in _CVD_MATRICES})
    minimum = {}
    for name, values in simulations.items():
        lab = _rgb_to_lab(values)
        minimum[name] = min(
            ciede2000(lab[first], lab[second])
            for first, second in itertools.combinations(range(len(lab)), 2)
        )
    linear = _srgb_to_linear(rgb)
    luminance = linear @ np.array([0.2126, 0.7152, 0.0722])
    return {
        "minimum_ciede2000": minimum,
        "grayscale_luminance": tuple(float(value) for value in luminance),
    }


__all__ = ["ciede2000", "palette_accessibility_report", "simulate_cvd"]
