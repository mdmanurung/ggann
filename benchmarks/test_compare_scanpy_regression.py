"""Tests for the fail-closed matched benchmark regression comparator."""

from __future__ import annotations

from copy import deepcopy

from benchmarks.compare_scanpy_regression import compare_documents


def _document() -> dict:
    result = {
        "case_id": "extended/csr/x/matrixplot",
        "parameters": {"workload": "matrixplot", "ggann_backend": "plotnine"},
        "selected_genes": ["g1", "g2"],
        "observed_shape": [100, 20],
        "input_bytes": {"total": 1000},
        "input_fingerprint_before": "same",
        "input_fingerprint_after": "same",
        "input_immutable": True,
        "comparability": {"status": "pass"},
        "stages": {
            stage: {"libraries": {"ggann": {"repeated": {"median_duration_seconds": value}}}}
            for stage, value in (
                ("preparation", 0.01),
                ("construction", 0.02),
                ("render", 0.45),
                ("end_to_end", 0.5),
            )
        },
        "isolated_memory": {
            "end_to_end": {
                "ggann": {
                    "summary": {
                        "peak_rss_delta_bytes": {"median": 10_000_000},
                        "retained_after_gc_bytes": {"median": 8_000_000},
                    }
                }
            }
        },
    }
    metadata = {
        field: value
        for field, value in {
            "python": "3.12",
            "platform": "test",
            "cpu_model": "cpu",
            "logical_cpus": 1,
            "packages": {"ggann": "0.1.0"},
            "thread_settings": {"OMP_NUM_THREADS": "1"},
            "figure": {"dpi": 80},
            "preset": "extended",
            "dataset": "synthetic",
            "formats": ["csr"],
            "workloads": ["matrixplot"],
            "sources": ["x"],
            "variants": ["base"],
            "shape_overrides": {},
            "repeats": 7,
            "ggann_backend": "plotnine",
            "seed": 7,
            "rss_interval_ms": 1.0,
            "isolated_memory_stages": ["end_to_end"],
            "isolated_memory_repeats": 7,
        }.items()
    }
    return {
        "schema_version": 1,
        "benchmark_kind": "ggann_scanpy_matched",
        "metadata": metadata,
        "results": [result],
    }


def test_comparator_passes_equal_or_faster_candidate():
    baseline = _document()
    candidate = deepcopy(baseline)
    candidate["results"][0]["stages"]["end_to_end"]["libraries"]["ggann"]["repeated"][
        "median_duration_seconds"
    ] = 0.4

    result = compare_documents(baseline, candidate)

    assert result["status"] == "pass"
    assert all(check["pass"] for check in result["checks"])


def test_comparator_fails_metric_above_five_percent():
    baseline = _document()
    candidate = deepcopy(baseline)
    candidate["results"][0]["isolated_memory"]["end_to_end"]["ggann"]["summary"][
        "peak_rss_delta_bytes"
    ]["median"] = 10_500_001

    result = compare_documents(baseline, candidate)

    assert result["status"] == "fail"
    assert any(not check["pass"] for check in result["checks"])


def test_comparator_fails_render_regression_above_five_percent():
    baseline = _document()
    candidate = deepcopy(baseline)
    candidate["results"][0]["stages"]["render"]["libraries"]["ggann"]["repeated"][
        "median_duration_seconds"
    ] = 0.48

    result = compare_documents(baseline, candidate)

    assert result["status"] == "fail"
    failed = [check for check in result["checks"] if not check["pass"]]
    assert [check["metric"] for check in failed] == ["render_median_seconds"]


def test_comparator_rejects_environment_mismatch():
    baseline = _document()
    candidate = deepcopy(baseline)
    candidate["metadata"]["packages"] = {"ggann": "0.2.0"}

    result = compare_documents(baseline, candidate)

    assert result["status"] == "not_comparable"
    assert result["issues"] == ["metadata.packages differs"]
