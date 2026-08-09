#!/usr/bin/env python3
"""Compare ggann timings and memory between matched Scanpy benchmark revisions."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any

_METADATA_FIELDS = (
    "python",
    "platform",
    "cpu_model",
    "logical_cpus",
    "packages",
    "thread_settings",
    "figure",
    "preset",
    "dataset",
    "formats",
    "workloads",
    "sources",
    "variants",
    "shape_overrides",
    "repeats",
    "ggann_backend",
    "seed",
    "rss_interval_ms",
    "isolated_memory_stages",
    "isolated_memory_repeats",
)
_CASE_FIELDS = (
    "parameters",
    "selected_genes",
    "observed_shape",
    "input_bytes",
    "input_fingerprint_before",
    "input_fingerprint_after",
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _case_map(document: dict[str, Any]) -> dict[str, dict[str, Any]]:
    return {result["case_id"]: result for result in document.get("results", [])}


def _comparability_issues(
    baseline: dict[str, Any], candidate: dict[str, Any]
) -> list[str]:
    issues = []
    for field in ("schema_version", "benchmark_kind"):
        if baseline.get(field) != candidate.get(field):
            issues.append(f"{field} differs")
    before_metadata = baseline.get("metadata", {})
    after_metadata = candidate.get("metadata", {})
    for field in _METADATA_FIELDS:
        if before_metadata.get(field) != after_metadata.get(field):
            issues.append(f"metadata.{field} differs")

    before_cases = _case_map(baseline)
    after_cases = _case_map(candidate)
    if set(before_cases) != set(after_cases):
        missing = sorted(set(before_cases) - set(after_cases))
        added = sorted(set(after_cases) - set(before_cases))
        if missing:
            issues.append("cases missing from candidate: " + ", ".join(missing))
        if added:
            issues.append("cases absent from baseline: " + ", ".join(added))
        return issues

    for case_id, before in before_cases.items():
        after = after_cases[case_id]
        for field in _CASE_FIELDS:
            if before.get(field) != after.get(field):
                issues.append(f"{case_id}.{field} differs")
        if not before.get("input_immutable") or not after.get("input_immutable"):
            issues.append(f"{case_id} mutated AnnData")
        if before.get("comparability", {}).get("status") != "pass":
            issues.append(f"{case_id} baseline payload is not comparable to Scanpy")
        if after.get("comparability", {}).get("status") != "pass":
            issues.append(f"{case_id} candidate payload is not comparable to Scanpy")
    return issues


def _ratio(candidate: float, baseline: float) -> float:
    if baseline:
        return candidate / baseline
    return 1.0 if not candidate else math.inf


def _timing_metric(result: dict[str, Any], stage: str) -> float:
    return float(
        result["stages"][stage]["libraries"]["ggann"]["repeated"][
            "median_duration_seconds"
        ]
    )


def _memory_metric(result: dict[str, Any], metric: str) -> float:
    return float(
        result["isolated_memory"]["end_to_end"]["ggann"]["summary"][metric]["median"]
    )


def compare_documents(
    baseline: dict[str, Any],
    candidate: dict[str, Any],
    *,
    max_regression_pct: float = 5.0,
) -> dict[str, Any]:
    """Evaluate the frozen-candidate timing and memory regression gate."""
    issues = _comparability_issues(baseline, candidate)
    if issues:
        return {
            "status": "not_comparable",
            "max_regression_pct": max_regression_pct,
            "issues": issues,
            "checks": [],
        }

    limit = 1 + max_regression_pct / 100
    before_cases = _case_map(baseline)
    after_cases = _case_map(candidate)
    checks = []
    for case_id in sorted(before_cases):
        workload = before_cases[case_id]["parameters"]["workload"]
        for stage in ("preparation", "end_to_end"):
            before = _timing_metric(before_cases[case_id], stage)
            after = _timing_metric(after_cases[case_id], stage)
            ratio = _ratio(after, before)
            checks.append(
                {
                    "case_id": case_id,
                    "workload": workload,
                    "metric": f"{stage}_median_seconds",
                    "kind": "timing",
                    "baseline": before,
                    "candidate": after,
                    "candidate_over_baseline": ratio,
                    "pass": ratio <= limit,
                }
            )
        for metric in ("peak_rss_delta_bytes", "retained_after_gc_bytes"):
            before = _memory_metric(before_cases[case_id], metric)
            after = _memory_metric(after_cases[case_id], metric)
            ratio = _ratio(after, before)
            checks.append(
                {
                    "case_id": case_id,
                    "workload": workload,
                    "metric": f"end_to_end_{metric}",
                    "kind": "memory",
                    "baseline": before,
                    "candidate": after,
                    "candidate_over_baseline": ratio,
                    "pass": ratio <= limit,
                }
            )
    return {
        "status": "pass" if all(check["pass"] for check in checks) else "fail",
        "max_regression_pct": max_regression_pct,
        "issues": [],
        "checks": checks,
    }


def _value(check: dict[str, Any], key: str) -> str:
    value = check[key]
    if check["kind"] == "timing":
        return f"{value * 1_000:.3f} ms"
    return f"{value / (1024 * 1024):.3f} MiB"


def _markdown(
    result: dict[str, Any],
    baseline_path: Path,
    candidate_path: Path,
    baseline: dict[str, Any],
    candidate: dict[str, Any],
) -> str:
    lines = [
        "# Frozen ggann candidate regression check",
        "",
        f"Verdict: **{result['status'].upper()}**",
        "",
        f"Maximum permitted regression: {result['max_regression_pct']:.1f}%.",
        f"Baseline: `{baseline_path}` (`{_sha256(baseline_path)}`).",
        f"Candidate: `{candidate_path}` (`{_sha256(candidate_path)}`).",
        "",
        "| Workload | Metric | Frozen | Candidate | Candidate / frozen | Result |",
        "|---|---|---:|---:|---:|:---:|",
    ]
    for check in result["checks"]:
        ratio = check["candidate_over_baseline"]
        ratio_text = "inf" if not math.isfinite(ratio) else f"{ratio:.3f}x"
        lines.append(
            f"| {check['workload']} | `{check['metric']}` | "
            f"{_value(check, 'baseline')} | {_value(check, 'candidate')} | "
            f"{ratio_text} | {'PASS' if check['pass'] else 'FAIL'} |"
        )
    if result["issues"]:
        lines.extend(["", "## Comparability failures", ""])
        lines.extend(f"- {issue}" for issue in result["issues"])
    lines.extend(
        [
            "",
            "## Source trees",
            "",
            f"- Frozen: `{baseline['metadata']['ggann_source']['python_tree_sha256']}`",
            f"- Candidate: `{candidate['metadata']['ggann_source']['python_tree_sha256']}`",
            "",
        ]
    )
    return "\n".join(lines)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("baseline", type=Path)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--max-regression-pct", type=float, default=5.0)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--json-output", type=Path)
    return parser


def main() -> int:
    args = _parser().parse_args()
    baseline = json.loads(args.baseline.read_text())
    candidate = json.loads(args.candidate.read_text())
    result = compare_documents(
        baseline,
        candidate,
        max_regression_pct=args.max_regression_pct,
    )
    report = _markdown(result, args.baseline, args.candidate, baseline, candidate)
    print(report, end="")
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(report)
    if args.json_output is not None:
        args.json_output.parent.mkdir(parents=True, exist_ok=True)
        args.json_output.write_text(json.dumps(result, indent=2) + "\n")
    return {"pass": 0, "fail": 1, "not_comparable": 2}[result["status"]]


if __name__ == "__main__":
    raise SystemExit(main())
