"""Focused tests for benchmark provenance and comparison guards."""

from __future__ import annotations

import unittest
from copy import deepcopy
from dataclasses import replace

import pandas as pd
from plotnine import aes, geom_point, ggplot, labs, scale_color_continuous

from benchmarks.check_invariance import _plot_contract, _provenance
from benchmarks.compare_results import _comparability_issues
from benchmarks.run_benchmarks import (
    CaseSpec,
    _frame_description,
    _ggann_source_metadata,
    _make_fixture,
    _selected_formats,
    _selected_workloads,
    _shape_variants,
    _stage_sizes,
    _workload_function,
)


def _document() -> dict:
    parameters = {
        "preset": "smoke",
        "matrix_format": "csr",
        "workload": "resolve_x",
        "n_obs": 10,
        "n_vars": 5,
        "n_genes": 2,
        "n_groups": 2,
        "density": 0.1,
        "embedding_dims": 2,
        "render_cells": 10,
        "seed": 7,
        "repeats": 2,
        "rss_interval_seconds": 0.001,
    }
    return {
        "schema_version": 1,
        "metadata": {
            "python": "3.x",
            "platform": "test-platform",
            "rss_backend": "proc_statm",
            "packages": {"ggann": "0.1.0", "numpy": "2.0"},
            "preset": "smoke",
            "repeats": 2,
            "seed": 7,
            "rss_interval_ms": 1.0,
        },
        "results": [
            {
                "case_id": "smoke/csr/resolve_x",
                "parameters": parameters,
                "selected_genes": ["gene_0", "gene_4"],
                "input_bytes": {"total": 100},
            }
        ],
    }


class ComparabilityTests(unittest.TestCase):
    def test_matching_legacy_documents_are_comparable(self) -> None:
        baseline = _document()
        candidate = deepcopy(baseline)
        candidate["metadata"]["packages"]["ggann"] = "0.2.0"
        self.assertEqual(_comparability_issues(baseline, candidate), [])

    def test_environment_mismatch_is_reported(self) -> None:
        baseline = _document()
        candidate = deepcopy(baseline)
        candidate["metadata"]["seed"] = 8
        candidate["metadata"]["packages"]["numpy"] = "2.1"
        issues = _comparability_issues(baseline, candidate)
        self.assertTrue(any("metadata.seed" in issue for issue in issues))
        self.assertTrue(any("metadata.packages" in issue for issue in issues))

    def test_fixture_mismatch_is_reported(self) -> None:
        baseline = _document()
        candidate = deepcopy(baseline)
        candidate["results"][0]["parameters"]["density"] = 0.2
        issues = _comparability_issues(baseline, candidate)
        self.assertIn("smoke/csr/resolve_x.parameters: fixture inputs differ", issues)

    def test_case_set_mismatch_is_reported(self) -> None:
        baseline = _document()
        candidate = deepcopy(baseline)
        candidate["results"] = []
        issues = _comparability_issues(baseline, candidate)
        self.assertEqual(issues, ["cases missing from candidate: smoke/csr/resolve_x"])


class ProvenanceTests(unittest.TestCase):
    def test_source_digest_is_stable(self) -> None:
        first = _ggann_source_metadata()
        second = _ggann_source_metadata()
        self.assertEqual(first, second)
        self.assertEqual(len(first["python_tree_sha256"]), 64)
        self.assertGreater(first["python_files"], 0)

    def test_artifact_provenance_records_source_and_environment(self) -> None:
        provenance = _provenance()
        self.assertEqual(len(provenance["ggann_source"]["python_tree_sha256"]), 64)
        self.assertIn("plotnine", provenance["packages"])
        self.assertEqual(
            set(provenance["thread_settings"]),
            {"OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"},
        )


class ArtifactContractTests(unittest.TestCase):
    def test_plot_contract_records_mapping_labels_and_scales(self) -> None:
        frame = pd.DataFrame({"x": [0.0, 1.0], "y": [1.0, 0.0]})
        plot = (
            ggplot(frame, aes("x", "y", color="x"))
            + geom_point()
            + scale_color_continuous()
            + labs(x="horizontal", color="signal")
        )

        contract = _plot_contract(plot)

        self.assertEqual(contract["mapping"], {"color": "x", "x": "x", "y": "y"})
        self.assertEqual(contract["labels"]["x"], "horizontal")
        self.assertEqual(contract["labels"]["color"], "signal")
        self.assertEqual(contract["scales"][0]["aesthetics"], ["color"])


class ScalingTests(unittest.TestCase):
    def test_each_shape_dimension_can_scale_in_isolation(self) -> None:
        expected = {
            "n_obs": [250, 500],
            "n_vars": [100, 200],
            "n_genes": [2, 4],
            "n_groups": [2, 4],
        }
        for field, values in expected.items():
            overrides = {name: None for name in expected}
            overrides[field] = values
            shapes = _shape_variants("smoke", "csr", overrides)
            self.assertEqual([shape[field] for shape in shapes], values)

    def test_multiple_scaling_dimensions_are_rejected(self) -> None:
        overrides = {
            "n_obs": [250, 500],
            "n_vars": [100, 200],
            "n_genes": None,
            "n_groups": None,
        }
        with self.assertRaisesRegex(ValueError, "Only one"):
            _shape_variants("smoke", "csr", overrides)

    def test_scaled_case_id_includes_changed_dimensions(self) -> None:
        spec = CaseSpec(
            preset="smoke",
            matrix_format="csr",
            workload="resolve_x",
            n_obs=500,
            n_vars=200,
            n_genes=4,
            n_groups=4,
            density=0.03,
            embedding_dims=10,
            render_cells=300,
            seed=7,
            repeats=1,
            rss_interval_seconds=0.001,
        )
        self.assertEqual(spec.case_id, "smoke/csr/resolve_x")
        scaled = replace(spec, n_obs=250, n_genes=2)
        self.assertEqual(scaled.case_id, "smoke[n_obs=250,n_genes=2]/csr/resolve_x")


class SelectionTests(unittest.TestCase):
    def test_empty_format_selection_is_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "format"):
            _selected_formats("")

    def test_empty_workload_selection_is_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "workload"):
            _selected_workloads("", include_rendering=False)


class StageSizeTests(unittest.TestCase):
    def test_core_stage_sizes_include_projection_and_prepared_data(self) -> None:
        spec = CaseSpec(
            preset="smoke",
            matrix_format="csr",
            workload="resolve_x",
            n_obs=10,
            n_vars=5,
            n_genes=2,
            n_groups=2,
            density=0.2,
            embedding_dims=2,
            render_cells=10,
            seed=7,
            repeats=1,
            rss_interval_seconds=0.001,
        )
        adata, genes, _ = _make_fixture(spec)
        output = _workload_function(adata, genes, spec.workload, 10)()
        description = _frame_description(output, fingerprint=True)
        stages = _stage_sizes(adata, genes, spec.workload, description)

        self.assertEqual(stages["projected_expression"]["shape"], [10, 2])
        self.assertGreater(stages["projected_expression"]["bytes"], 0)
        self.assertEqual(stages["projected_obs"]["shape"], [10, 2])
        self.assertEqual(stages["prepared_data"]["shape"], [10, 4])
        self.assertEqual(stages["prepared_data"]["bytes"], description["bytes"])


if __name__ == "__main__":
    unittest.main()
