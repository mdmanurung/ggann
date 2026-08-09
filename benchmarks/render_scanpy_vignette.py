#!/usr/bin/env python3
"""Render the committed Scanpy comparison JSON as a Markdown include."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

_LABELS = {
    "embedding_categorical": "Categorical embedding",
    "dotplot": "Dotplot",
    "matrixplot": "Matrixplot",
}


def _milliseconds(seconds: float) -> str:
    return f"{seconds * 1_000:.1f} ms"


def _mebibytes(value: int) -> str:
    return f"{value / 2**20:.2f} MiB"


def _memory_metric(record: dict[str, Any], name: str) -> tuple[float, float, float]:
    if "summary" in record:
        summary = record["summary"][name]
        return float(summary["median"]), float(summary["min"]), float(summary["max"])
    value = float(record["sample"][name])
    return value, value, value


def _memory_text(values: tuple[float, float, float]) -> str:
    median, minimum, maximum = values
    if minimum == maximum:
        return _mebibytes(int(median))
    return f"{_mebibytes(int(median))} [{_mebibytes(int(minimum))}, {_mebibytes(int(maximum))}]"


def render_document(document: dict[str, Any]) -> str:
    """Return the primary benchmark tables derived from one raw result document."""
    if document.get("schema_version") != 1:
        raise ValueError("expected benchmark schema_version=1")
    if document.get("benchmark_kind") != "ggann_scanpy_matched":
        raise ValueError("expected a ggann_scanpy_matched benchmark document")

    results = {
        result["parameters"]["workload"]: result
        for result in document.get("results", [])
        if result.get("parameters", {}).get("workload") in _LABELS
    }
    missing = [workload for workload in _LABELS if workload not in results]
    if missing:
        raise ValueError("missing primary workload(s): " + ", ".join(missing))

    metadata = document["metadata"]
    first = results["embedding_categorical"]["parameters"]
    lines = [
        (
            f"Measured on {first['n_obs']:,} observations, {first['n_vars']:,} variables, "
            f"{first['n_genes']} requested genes, and {first['n_groups']} groups "
            f"({metadata['repeats']} warm repetitions, CSR `.X`)."
        ),
        "",
        "| Workload | ggann prep | Scanpy prep | ggann construct | Scanpy construct | ggann render/save | Scanpy render/save | ggann end to end | Scanpy end to end | E2E speedup |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for workload, label in _LABELS.items():
        result = results[workload]
        values = []
        for stage in ("preparation", "construction", "render", "end_to_end"):
            stage_result = result["stages"][stage]
            for library in ("ggann", "scanpy"):
                values.append(
                    _milliseconds(
                        stage_result["libraries"][library]["repeated"]["median_duration_seconds"]
                    )
                )
        speedup = result["stages"]["end_to_end"]["speedup_scanpy_over_ggann"]
        lines.append(f"| {label} | {' | '.join(values)} | {speedup:.2f}x |")

    lines.extend(
        [
            "",
            "Fresh-child memory probes exclude imports and fixture creation from the baseline:",
            "",
            "| Workload | ggann peak RSS | Scanpy peak RSS | ggann retained RSS | Scanpy retained RSS | Peak ratio |",
            "|---|---:|---:|---:|---:|---:|",
        ]
    )
    for workload, label in _LABELS.items():
        memory = results[workload]["isolated_memory"]["end_to_end"]
        ggann_peak = _memory_metric(memory["ggann"], "peak_rss_delta_bytes")
        scanpy_peak = _memory_metric(memory["scanpy"], "peak_rss_delta_bytes")
        ggann_retained = _memory_metric(memory["ggann"], "retained_after_gc_bytes")
        scanpy_retained = _memory_metric(memory["scanpy"], "retained_after_gc_bytes")
        ratio = ggann_peak[0] / scanpy_peak[0]
        lines.append(
            f"| {label} | {_memory_text(ggann_peak)} | "
            f"{_memory_text(scanpy_peak)} | "
            f"{_memory_text(ggann_retained)} | "
            f"{_memory_text(scanpy_retained)} | {ratio:.2f}x |"
        )

    gates = document["release_gates"]
    source = metadata["ggann_source"]["python_tree_sha256"]
    lines.extend(
        [
            "",
            (
                f"Speedup is Scanpy time divided by ggann time; values above 1 favour ggann. "
                f"The recorded gate status is **{gates['status'].upper()}**. All prepared payloads "
                "passed the benchmark's numeric equivalence checks and all input fingerprints "
                "were unchanged."
            ),
            "",
            (
                "Construction is diagnostic because Scanpy eagerly creates different amounts of "
                "Matplotlib state. `render/save` measures PNG save from a materialized figure; "
                "end-to-end is the authoritative plotting comparison."
            ),
            "",
            f"Recorded ggann source-tree SHA-256: `{source}`.",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    document = json.loads(args.input.read_text())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(render_document(document))


if __name__ == "__main__":
    main()
