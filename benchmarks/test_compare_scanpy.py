"""Focused tests for the matched ggann-versus-Scanpy benchmark runner."""

from __future__ import annotations

import unittest

import pandas as pd

from benchmarks.compare_scanpy import (
    ComparisonSpec,
    _adata_digest,
    _case_supported,
    _compare_prepared,
    _evaluate_release_gates,
    _execute_case,
    _ggann_preparation,
    _scanpy_preparation,
    _selected_sources,
    _selected_workloads,
)
from benchmarks.run_benchmarks import _make_fixture


def _spec(workload: str = "dotplot") -> ComparisonSpec:
    return ComparisonSpec(
        preset="smoke",
        matrix_format="csr",
        workload=workload,
        source="x",
        n_obs=24,
        n_vars=12,
        n_genes=3,
        n_groups=3,
        density=0.2,
        embedding_dims=4,
        render_cells=24,
        seed=7,
        repeats=1,
        rss_interval_seconds=0.001,
    )


class SelectionTests(unittest.TestCase):
    def test_aliases_expand_to_documented_workloads_and_sources(self) -> None:
        self.assertEqual(
            _selected_workloads("primary"),
            ["embedding_categorical", "dotplot", "matrixplot"],
        )
        self.assertEqual(_selected_sources("all"), ["x", "layer", "raw"])

    def test_raw_highest_expression_is_explicitly_unsupported(self) -> None:
        supported, reason = _case_supported("highest_expr_genes", "raw")
        self.assertFalse(supported)
        self.assertIn("no use_raw", reason)


class ComparabilityTests(unittest.TestCase):
    def test_numeric_tolerance_and_categorical_color_are_supported(self) -> None:
        left = pd.DataFrame(
            {
                "obs_name": ["a", "b"],
                "x": [0.0, 1.0],
                "y": [1.0, 0.0],
                "color": pd.Categorical(["T", "B"], categories=["T", "B"]),
            }
        )
        right = pd.DataFrame(
            {
                "obs_name": ["a", "b"],
                "x": [0.0, 1.0 + 1e-8],
                "y": [1.0, 0.0],
                "color": pd.Categorical(["T", "B"], categories=["B", "T"]),
            }
        )
        result = _compare_prepared(left, right, "embedding_categorical")
        self.assertEqual(result["status"], "pass")

    def test_dotplot_preparation_matches_on_sparse_fixture(self) -> None:
        spec = _spec("dotplot")
        adata, genes, _ = _make_fixture(spec.fixture_spec())
        left = _ggann_preparation(adata, genes, spec)
        right = _scanpy_preparation(adata, genes, spec)
        result = _compare_prepared(left, right, spec.workload)
        self.assertEqual(result["status"], "pass", result["issues"])


class InputOwnershipTests(unittest.TestCase):
    def test_digest_changes_after_expression_mutation(self) -> None:
        spec = _spec()
        adata, _, _ = _make_fixture(spec.fixture_spec())
        before = _adata_digest(adata)
        adata.X.data[0] += 1
        self.assertNotEqual(before, _adata_digest(adata))


class RunnerSmokeTests(unittest.TestCase):
    def test_small_embedding_case_records_all_stages(self) -> None:
        result = _execute_case(_spec("embedding_categorical"))
        self.assertEqual(result["comparability"]["status"], "pass")
        self.assertTrue(result["input_immutable"])
        self.assertEqual(
            set(result["stages"]),
            {"preparation", "construction", "render", "end_to_end"},
        )
        self.assertEqual(result["stages"]["construction"]["claim_status"], "diagnostic_only")
        for stage in result["stages"].values():
            for library in ("ggann", "scanpy"):
                self.assertEqual(len(stage["libraries"][library]["repeated"]["samples"]), 1)

    def test_release_gates_fail_closed_without_large_primary_cases(self) -> None:
        gates = _evaluate_release_gates([])
        self.assertEqual(gates["status"], "not_evaluated")
        self.assertEqual(
            gates["baseline_regression_gate"],
            "not_evaluated_by_this_cross-library_runner",
        )


if __name__ == "__main__":
    unittest.main()
