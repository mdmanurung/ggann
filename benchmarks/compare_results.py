#!/usr/bin/env python3
"""Compare two JSON outputs from run_benchmarks.py as a Markdown table."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys
from typing import Any


def _load(path: Path) -> dict[str, Any]:
    document = json.loads(path.read_text())
    if document.get("schema_version") != 1 or not isinstance(
        document.get("results"), list
    ):
        raise ValueError(f"{path} is not a ggann benchmark schema-version 1 document")
    return document


def _index(document: dict[str, Any]) -> dict[str, dict[str, Any]]:
    return {result["case_id"]: result for result in document["results"]}


def _comparability_issues(
    baseline: dict[str, Any], candidate: dict[str, Any]
) -> list[str]:
    """Return reasons two benchmark documents are not directly comparable."""
    issues: list[str] = []
    before_metadata = baseline.get("metadata", {})
    after_metadata = candidate.get("metadata", {})

    for field in (
        "python",
        "platform",
        "rss_backend",
        "preset",
        "repeats",
        "seed",
        "rss_interval_ms",
    ):
        before = before_metadata.get(field)
        after = after_metadata.get(field)
        if before != after:
            issues.append(f"metadata.{field}: {before!r} != {after!r}")

    before_packages = dict(before_metadata.get("packages", {}))
    after_packages = dict(after_metadata.get("packages", {}))
    # ggann is the subject of the comparison; its version may legitimately change.
    before_packages.pop("ggann", None)
    after_packages.pop("ggann", None)
    if before_packages != after_packages:
        issues.append("metadata.packages: dependency versions differ")

    before_threads = before_metadata.get("thread_settings")
    after_threads = after_metadata.get("thread_settings")
    # Schema-version 1 documents created before thread provenance was added do
    # not contain this field. Preserve comparisons between those documents.
    if (
        before_threads is not None
        and after_threads is not None
        and before_threads != after_threads
    ):
        issues.append(
            f"metadata.thread_settings: {before_threads!r} != {after_threads!r}"
        )

    before_by_id = _index(baseline)
    after_by_id = _index(candidate)
    before_ids = set(before_by_id)
    after_ids = set(after_by_id)
    missing = sorted(before_ids - after_ids)
    added = sorted(after_ids - before_ids)
    if missing:
        issues.append("cases missing from candidate: " + ", ".join(missing))
    if added:
        issues.append("cases only in candidate: " + ", ".join(added))

    for case_id in sorted(before_ids & after_ids):
        before = before_by_id[case_id]
        after = after_by_id[case_id]
        for field in ("parameters", "selected_genes", "input_bytes"):
            if before.get(field) != after.get(field):
                issues.append(f"{case_id}.{field}: fixture inputs differ")
    return issues


def _percent_change(before: float, after: float) -> float:
    if before == 0:
        return 0.0 if after == 0 else math.inf
    return (after - before) / before * 100


def _percent_text(before: float, after: float) -> str:
    change = _percent_change(before, after)
    if math.isinf(change):
        return "+inf%"
    return f"{change:+.1f}%"


def _seconds(value: float) -> str:
    if value < 0.001:
        return f"{value * 1_000_000:.0f} us"
    if value < 1:
        return f"{value * 1_000:.1f} ms"
    return f"{value:.3f} s"


def _bytes(value: int) -> str:
    units = ("B", "KiB", "MiB", "GiB")
    amount = float(value)
    for unit in units:
        if amount < 1024 or unit == units[-1]:
            return f"{amount:.1f} {unit}"
        amount /= 1024
    raise AssertionError("unreachable")


def _escape(value: str) -> str:
    return value.replace("|", "\\|")


def _table(
    baseline: dict[str, Any], candidate: dict[str, Any]
) -> tuple[str, list[dict[str, Any]]]:
    before_by_id = _index(baseline)
    after_by_id = _index(candidate)
    shared = sorted(before_by_id.keys() & after_by_id.keys())
    rows = []
    lines = [
        "| Case | Cold before | Cold after | Cold change | Repeat before | Repeat after | Repeat change | Peak RSS before | Peak RSS after | Peak change | Retained before | Retained after | Output |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|:---:|",
    ]
    for case_id in shared:
        before = before_by_id[case_id]
        after = after_by_id[case_id]
        cold_before = before["cold"]["duration_seconds"]
        cold_after = after["cold"]["duration_seconds"]
        repeat_before = before["repeated"]["median_duration_seconds"]
        repeat_after = after["repeated"]["median_duration_seconds"]
        peak_before = max(
            before["cold"]["peak_rss_delta_bytes"],
            before["repeated"]["max_peak_rss_delta_bytes"],
        )
        peak_after = max(
            after["cold"]["peak_rss_delta_bytes"],
            after["repeated"]["max_peak_rss_delta_bytes"],
        )
        retained_before = max(
            before["cold"]["retained_after_gc_bytes"],
            before["repeated"]["max_retained_after_gc_bytes"],
        )
        retained_after = max(
            after["cold"]["retained_after_gc_bytes"],
            after["repeated"]["max_retained_after_gc_bytes"],
        )
        output_equal = before["output"]["fingerprint"] == after["output"]["fingerprint"]
        row = {
            "case_id": case_id,
            "cold_change_pct": _percent_change(cold_before, cold_after),
            "repeat_change_pct": _percent_change(repeat_before, repeat_after),
            "peak_change_pct": _percent_change(peak_before, peak_after),
            "retained_change_pct": _percent_change(retained_before, retained_after),
            "output_equal": output_equal,
        }
        rows.append(row)
        lines.append(
            "| "
            + " | ".join(
                (
                    _escape(case_id),
                    _seconds(cold_before),
                    _seconds(cold_after),
                    _percent_text(cold_before, cold_after),
                    _seconds(repeat_before),
                    _seconds(repeat_after),
                    _percent_text(repeat_before, repeat_after),
                    _bytes(peak_before),
                    _bytes(peak_after),
                    _percent_text(peak_before, peak_after),
                    _bytes(retained_before),
                    _bytes(retained_after),
                    "same" if output_equal else "CHANGED",
                )
            )
            + " |"
        )

    missing_after = sorted(before_by_id.keys() - after_by_id.keys())
    new_after = sorted(after_by_id.keys() - before_by_id.keys())
    if missing_after:
        lines.extend(("", "Missing from candidate: " + ", ".join(missing_after)))
    if new_after:
        lines.extend(("", "Only in candidate: " + ", ".join(new_after)))
    return "\n".join(lines) + "\n", rows


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("baseline", type=Path)
    parser.add_argument("candidate", type=Path)
    parser.add_argument(
        "--output", type=Path, help="Also write the Markdown table here."
    )
    parser.add_argument(
        "--fail-regression-pct",
        type=float,
        help="Fail when repeated time or peak RSS regresses by more than this percentage.",
    )
    parser.add_argument(
        "--fail-on-output-change",
        action="store_true",
        help="Fail when the prepared output fingerprint changes.",
    )
    parser.add_argument(
        "--allow-incomparable",
        action="store_true",
        help=(
            "Compare mismatched environments or fixtures for exploratory use; "
            "the report will be marked as incomparable."
        ),
    )
    return parser


def main() -> int:
    args = _build_parser().parse_args()
    baseline = _load(args.baseline)
    candidate = _load(args.candidate)
    issues = _comparability_issues(baseline, candidate)
    if issues and not args.allow_incomparable:
        print(
            "Benchmark inputs are not comparable:\n"
            + "\n".join(f"- {issue}" for issue in issues),
            file=sys.stderr,
        )
        print(
            "Use --allow-incomparable only for an explicitly exploratory report.",
            file=sys.stderr,
        )
        return 2
    table, rows = _table(baseline, candidate)

    before_label = baseline.get("metadata", {}).get("label", "baseline")
    after_label = candidate.get("metadata", {}).get("label", "candidate")
    heading = f"# ggann benchmark comparison: {before_label} vs {after_label}\n\n"
    warning = ""
    if issues:
        warning = (
            "> **Warning:** benchmark inputs are not comparable.\n>\n"
            + "\n".join(f"> - {issue}" for issue in issues)
            + "\n\n"
        )
    report = heading + warning + table
    print(report, end="")
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(report)

    failed = False
    if args.fail_regression_pct is not None:
        threshold = args.fail_regression_pct
        failed = any(
            row["repeat_change_pct"] > threshold or row["peak_change_pct"] > threshold
            for row in rows
        )
    if args.fail_on_output_change:
        failed = failed or any(not row["output_equal"] for row in rows)
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
