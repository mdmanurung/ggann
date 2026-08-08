"""Every plotnine-native utility helper must have a layer-by-layer grammar twin.

The twins live in ``examples/grammar_equivalents.py``. These tests build each one
so the documented grammar cannot silently drift from the helpers, and assert the
mapping covers every plotnine-native utility (the escape hatches -- plot_clustermap
and plot_upset -- are not plotnine and are intentionally excluded).
"""

from __future__ import annotations

import importlib.util
import os

import numpy as np
import plotnine as p9
import pytest

_GE_PATH = os.path.join(
    os.path.dirname(__file__), "..", "examples", "grammar_equivalents.py"
)
_spec = importlib.util.spec_from_file_location("grammar_equivalents", _GE_PATH)
ge = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ge)


@pytest.mark.parametrize("name,builder", list(ge.EQUIVALENTS.items()))
def test_grammar_equivalent_builds(adata, name, builder):
    plot = builder(adata)
    assert isinstance(plot, p9.ggplot)
    plot._build()  # force layout so any layer mismatch surfaces here


@pytest.mark.parametrize("name", list(ge.HELPERS))
def test_comparison_helper_builds(adata, name):
    # each grammar twin is paired with the one-line helper it reproduces
    plot = ge.HELPERS[name](adata)
    assert isinstance(plot, p9.ggplot)
    plot._build()


def test_every_plotnine_utility_has_a_grammar_equivalent():
    expected = {
        "plot_embedding",
        "plot_features",
        "plot_density",
        "plot_dotplot",
        "plot_matrixplot",
        "plot_violin",
        "plot_box",
        "plot_expression_bar",
        "plot_expression_line",
        "plot_proportions",
        "plot_correlation",
    }
    missing = expected - set(ge.EQUIVALENTS)
    assert not missing, f"utilities lacking a grammar equivalent: {sorted(missing)}"
    # every grammar twin must have a paired helper for the side-by-side comparison
    assert set(ge.EQUIVALENTS) == set(ge.HELPERS)


@pytest.mark.parametrize("has_raw", [True, False])
def test_correlation_reference_matches_helper_data(adata, has_raw):
    selected = adata if has_raw else adata.copy()
    if not has_raw:
        selected.raw = None

    helper = ge.ag.plot_correlation(selected, ge.GROUP, cluster=False).data.copy()
    reference = ge.correlation_grammar(selected).data.copy()

    def keyed(frame):
        frame["row"] = frame["row"].astype(str)
        frame["col"] = frame["col"].astype(str)
        return frame.set_index(["row", "col"])["corr"].sort_index()

    helper_values = keyed(helper)
    reference_values = keyed(reference)
    assert helper_values.index.equals(reference_values.index)
    np.testing.assert_allclose(helper_values, reference_values, rtol=1e-12, atol=1e-12)
