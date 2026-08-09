"""Accessibility gates for the reviewed publication colour vocabulary."""

from __future__ import annotations

import numpy as np
import pytest
from matplotlib.colors import to_rgb

import ggann as ag
from ggann._palette_qa import ciede2000, palette_accessibility_report


def _relative_luminance(colours) -> np.ndarray:
    rgb = np.array([to_rgb(colour) for colour in colours])
    linear = np.where(rgb <= 0.04045, rgb / 12.92, ((rgb + 0.055) / 1.055) ** 2.4)
    return linear @ np.array([0.2126, 0.7152, 0.0722])


def test_ciede2000_reference_pair():
    first = np.array([50.0, 2.6772, -79.7751])
    second = np.array([50.0, 0.0, -82.7485])
    assert ciede2000(first, second) == pytest.approx(2.0425, abs=1e-4)


def test_core_eight_palette_separation_normal_and_cvd():
    report = palette_accessibility_report(ag.publication_palette("qualitative", 8))
    separation = report["minimum_ciede2000"]
    assert separation["normal"] >= 10
    assert separation["protanopia"] >= 5
    assert separation["deuteranopia"] >= 5
    assert separation["tritanopia"] >= 5
    assert len(report["grayscale_luminance"]) == 8


def test_sequential_and_diverging_luminance_are_monotonic():
    sequential = _relative_luminance(ag.publication_palette("sequential", 33))
    assert np.all(np.diff(sequential) >= -1e-10)

    diverging = _relative_luminance(ag.publication_palette("diverging", 33))
    midpoint = len(diverging) // 2
    assert np.all(np.diff(diverging[: midpoint + 1]) >= -1e-10)
    assert np.all(np.diff(diverging[midpoint:]) <= 1e-10)
