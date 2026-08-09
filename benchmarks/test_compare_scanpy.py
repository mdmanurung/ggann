"""Focused tests for the matched ggann-versus-Scanpy benchmark runner."""

from __future__ import annotations

import unittest
import warnings
from dataclasses import replace

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import sparse

from benchmarks.compare_scanpy import (
    ComparisonSpec,
    _adata_digest,
    _case_supported,
    _compare_prepared,
    _evaluate_release_gates,
    _execute_case,
    _ggann_highest_frame,
    _ggann_native_preparation,
    _ggann_preparation,
    _scanpy_highest_frame,
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

    def test_matplotlib_backend_uses_compact_native_primary_payloads(self) -> None:
        for workload in ("embedding_categorical", "dotplot", "matrixplot"):
            with self.subTest(workload=workload):
                spec = replace(_spec(workload), ggann_backend="matplotlib")
                adata, genes, _ = _make_fixture(spec.fixture_spec())
                canonical = _ggann_preparation(adata, genes, spec)
                scanpy = _scanpy_preparation(adata, genes, spec)
                native = _ggann_native_preparation(adata, genes, spec)

                result = _compare_prepared(canonical, scanpy, workload)
                self.assertEqual(result["status"], "pass", result["issues"])
                self.assertIn("ggann_backend=matplotlib", spec.case_id)
                if workload == "embedding_categorical":
                    self.assertEqual(native.shape[1], 3)
                elif workload == "dotplot":
                    self.assertEqual(set(native), {"mean_expression", "fraction"})
                else:
                    self.assertIsInstance(native, pd.DataFrame)

    def test_highest_expression_zero_total_rows_match_scanpy(self) -> None:
        values = np.array(
            [
                [0.0, 0.0, 0.0],
                [1.0, 3.0, 0.0],
                [2.0, 0.0, 2.0],
            ]
        )
        formats = {
            "dense": lambda matrix: matrix,
            "csr": sparse.csr_matrix,
            "csc": sparse.csc_matrix,
        }
        for name, convert in formats.items():
            with self.subTest(matrix_format=name):
                adata = AnnData(
                    convert(values),
                    var=pd.DataFrame(index=["g1", "g2", "g3"]),
                )
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", UserWarning)
                    left = _ggann_highest_frame(adata, 3, "x")
                    right = _scanpy_highest_frame(adata, 3, "x")
                result = _compare_prepared(left, right, "highest_expr_genes")
                self.assertEqual(result["status"], "pass", result["issues"])

    def test_highest_expression_nan_difference_remains_visible(self) -> None:
        values = np.array(
            [
                [1.0, np.nan, 3.0],
                [2.0, 2.0, 0.0],
            ]
        )
        adata = AnnData(values, var=pd.DataFrame(index=["g1", "g2", "g3"]))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            left = _ggann_highest_frame(adata, 3, "x")
            right = _scanpy_highest_frame(adata, 3, "x")

        result = _compare_prepared(left, right, "highest_expr_genes")

        self.assertEqual(result["status"], "fail")
        self.assertIn("prepared numeric column 'percent' differs", result["issues"])


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

    def test_small_direct_backend_case_records_backend_and_passes(self) -> None:
        spec = replace(_spec("matrixplot"), ggann_backend="matplotlib")
        result = _execute_case(spec)
        self.assertEqual(result["comparability"]["status"], "pass")
        self.assertEqual(result["parameters"]["ggann_backend"], "matplotlib")
        self.assertTrue(result["input_immutable"])

    def test_release_gates_fail_closed_without_large_primary_cases(self) -> None:
        gates = _evaluate_release_gates([])
        self.assertEqual(gates["status"], "not_evaluated")
        self.assertEqual(
            gates["baseline_regression_gate"],
            "not_evaluated_by_this_cross-library_runner",
        )


if __name__ == "__main__":
    unittest.main()
