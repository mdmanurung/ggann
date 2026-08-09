#!/usr/bin/env python3
"""Profile concrete ggann preparation and plotnine rendering hot spots.

This diagnostic complements :mod:`benchmarks.compare_scanpy`: the comparison
runner supplies claim-quality timings, while this runner records call profiles,
Python allocations, RSS deltas, DataFrame sizes, and annplyr extraction calls
for the same deterministic fixture.  Timings collected under cProfile or
tracemalloc are deliberately labelled diagnostic and are never used as speed
claims.
"""

from __future__ import annotations

import argparse
import cProfile
import gc
import hashlib
import json
import os
import pstats
import statistics
import tempfile
import time
import tracemalloc
from contextlib import contextmanager
from dataclasses import asdict
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, Callable, Iterator

for _name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_name, "1")
os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "ggann-benchmark-mpl"))

try:
    from benchmarks.compare_scanpy import (
        ComparisonSpec,
        _adata_digest,
        _add_fixture_metadata,
        _close_output,
        _construct,
        _end_to_end_function,
        _ggann_native_preparation,
    )
    from benchmarks.run_benchmarks import (
        _ggann_source_metadata,
        _make_fixture,
        _rss_bytes,
        _RSSSampler,
        _shape_for,
    )
except ModuleNotFoundError:  # direct ``python benchmarks/profile_hotspots.py``
    from compare_scanpy import (  # type: ignore[no-redef]
        ComparisonSpec,
        _adata_digest,
        _add_fixture_metadata,
        _close_output,
        _construct,
        _end_to_end_function,
        _ggann_native_preparation,
    )
    from run_benchmarks import (  # type: ignore[no-redef]
        _ggann_source_metadata,
        _make_fixture,
        _rss_bytes,
        _RSSSampler,
        _shape_for,
    )


_WORKLOADS = ("embedding_categorical", "dotplot", "matrixplot")
_STAGES = ("preparation", "construction", "draw", "save", "end_to_end")


def _jsonable_selector(value: Any) -> Any:
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    if isinstance(value, dict):
        return {str(key): _jsonable_selector(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, range)):
        return [_jsonable_selector(item) for item in value]
    return repr(value)


def _frame_description(value: Any) -> dict[str, Any]:
    import pandas as pd

    if not isinstance(value, pd.DataFrame):
        description = {"kind": f"{type(value).__module__}.{type(value).__name__}"}
        plot_data = getattr(value, "data", None)
        if isinstance(plot_data, pd.DataFrame):
            description["plot_data"] = _frame_description(plot_data)
        return description
    sparse_columns = sum(isinstance(dtype, pd.SparseDtype) for dtype in value.dtypes)
    return {
        "kind": "pandas.DataFrame",
        "shape": [int(value.shape[0]), int(value.shape[1])],
        "bytes": int(value.memory_usage(index=True, deep=True).sum()),
        "sparse_columns": sparse_columns,
        "dtypes": [str(dtype) for dtype in value.dtypes],
    }


@contextmanager
def _record_annplyr_calls(adata: Any) -> Iterator[list[dict[str, Any]]]:
    """Record public ``to_df`` calls without bypassing annplyr semantics."""
    accessor_type = type(adata.ap)
    original = accessor_type.to_df
    calls: list[dict[str, Any]] = []

    def wrapped(self: Any, *args: Any, **kwargs: Any) -> Any:
        started = time.perf_counter_ns()
        result = original(self, *args, **kwargs)
        calls.append(
            {
                "duration_seconds": (time.perf_counter_ns() - started) / 1_000_000_000,
                "arguments": {
                    "positional": _jsonable_selector(args),
                    **{key: _jsonable_selector(value) for key, value in kwargs.items()},
                },
                "result": _frame_description(result),
            }
        )
        return result

    accessor_type.to_df = wrapped
    try:
        yield calls
    finally:
        accessor_type.to_df = original


def _profile_rows(profile: cProfile.Profile, *, limit: int = 80) -> list[dict[str, Any]]:
    stats = pstats.Stats(profile)
    rows = []
    for (filename, line, name), (
        primitive,
        calls,
        total,
        cumulative,
        _callers,
    ) in stats.stats.items():
        rows.append(
            {
                "file": filename,
                "line": line,
                "function": name,
                "primitive_calls": primitive,
                "calls": calls,
                "self_seconds": total,
                "cumulative_seconds": cumulative,
            }
        )
    rows.sort(key=lambda row: (row["cumulative_seconds"], row["self_seconds"]), reverse=True)
    return rows[:limit]


def _allocation_rows(snapshot: tracemalloc.Snapshot, *, limit: int = 40) -> list[dict[str, Any]]:
    rows = []
    for statistic in snapshot.statistics("lineno")[:limit]:
        frame = statistic.traceback[0]
        rows.append(
            {
                "file": frame.filename,
                "line": frame.lineno,
                "bytes": statistic.size,
                "blocks": statistic.count,
            }
        )
    return rows


def _stage_call(
    stage: str,
    adata: Any,
    genes: list[str],
    spec: ComparisonSpec,
    output_path: Path,
) -> tuple[Callable[[], Any], Callable[[Any], None]]:
    from matplotlib import pyplot as plt

    if stage == "preparation":
        return lambda: _ggann_native_preparation(adata, genes, spec), lambda _output: None
    if stage == "construction":
        return lambda: _construct("ggann", adata, genes, spec), _close_output
    if stage == "draw":
        plot = _construct("ggann", adata, genes, spec)

        def cleanup(figure: Any) -> None:
            plt.close(figure)
            _close_output(plot)

        return lambda: plot.draw(show=False), cleanup
    if stage == "save":
        plot = _construct("ggann", adata, genes, spec)
        figure = plot.draw(show=False)
        figure.set_size_inches(6.0, 4.5, forward=True)

        def save() -> dict[str, int]:
            figure.savefig(output_path, format="png", dpi=80)
            return {"png_bytes": output_path.stat().st_size, "axes": len(figure.axes)}

        def cleanup(_result: Any) -> None:
            plt.close(figure)
            _close_output(plot)

        return save, cleanup
    if stage == "end_to_end":
        return _end_to_end_function("ggann", adata, genes, spec, output_path), lambda _output: None
    raise ValueError(f"Unknown stage: {stage}")


def _measure_stage(
    stage: str,
    adata: Any,
    genes: list[str],
    spec: ComparisonSpec,
    output_dir: Path,
    profile_dir: Path,
) -> dict[str, Any]:
    samples = []
    for repeat in range(spec.repeats):
        function, cleanup = _stage_call(
            stage, adata, genes, spec, output_dir / f"{spec.workload}-{stage}-{repeat}.png"
        )
        gc.collect()
        baseline = _rss_bytes()
        with _record_annplyr_calls(adata) as extraction_calls:
            with _RSSSampler(spec.rss_interval_seconds) as sampler:
                started = time.perf_counter_ns()
                result = function()
                duration = (time.perf_counter_ns() - started) / 1_000_000_000
        with_output = _rss_bytes()
        description = _frame_description(result)
        cleanup(result)
        del result
        gc.collect()
        after_gc = _rss_bytes()
        samples.append(
            {
                "duration_seconds": duration,
                "baseline_rss_bytes": baseline,
                "peak_rss_delta_bytes": max(0, sampler.peak_bytes - baseline),
                "retained_with_output_bytes": max(0, with_output - baseline),
                "retained_after_gc_bytes": max(0, after_gc - baseline),
                "output": description,
                "annplyr_calls": extraction_calls,
            }
        )

    function, cleanup = _stage_call(
        stage, adata, genes, spec, output_dir / f"{spec.workload}-{stage}-profile.png"
    )
    profile = cProfile.Profile()
    with _record_annplyr_calls(adata) as profile_extractions:
        profiled_result = profile.runcall(function)
    profile_path = profile_dir / f"{spec.workload}-{stage}.prof"
    profile.dump_stats(profile_path)
    cleanup(profiled_result)
    del profiled_result
    gc.collect()

    function, cleanup = _stage_call(
        stage, adata, genes, spec, output_dir / f"{spec.workload}-{stage}-allocations.png"
    )
    tracemalloc.start()
    allocation_result = function()
    current_bytes, peak_bytes = tracemalloc.get_traced_memory()
    snapshot = tracemalloc.take_snapshot()
    tracemalloc.stop()
    cleanup(allocation_result)
    del allocation_result
    gc.collect()

    durations = [sample["duration_seconds"] for sample in samples]
    return {
        "stage": stage,
        "timing_scope": "unprofiled warm in-process diagnostic",
        "samples": samples,
        "median_duration_seconds": statistics.median(durations),
        "min_duration_seconds": min(durations),
        "max_duration_seconds": max(durations),
        "profile": {
            "path": str(profile_path),
            "sha256": hashlib.sha256(profile_path.read_bytes()).hexdigest(),
            "annplyr_calls": profile_extractions,
            "top_functions": _profile_rows(profile),
        },
        "allocations": {
            "current_bytes": current_bytes,
            "peak_bytes": peak_bytes,
            "top_lines": _allocation_rows(snapshot),
        },
    }


def _metadata(args: argparse.Namespace) -> dict[str, Any]:
    import importlib.metadata
    import platform
    import sys

    packages = {}
    for package in (
        "ggann",
        "scanpy",
        "annplyr",
        "anndata",
        "numpy",
        "pandas",
        "scipy",
        "plotnine",
        "plotnine-extra",
        "matplotlib",
    ):
        try:
            packages[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            packages[package] = None
    cpuinfo = Path("/proc/cpuinfo")
    cpu_model = platform.processor()
    if cpuinfo.exists():
        cpu_model = next(
            (
                line.partition(":")[2].strip()
                for line in cpuinfo.read_text().splitlines()
                if line.startswith("model name")
            ),
            cpu_model,
        )
    return {
        "created_at": datetime.now(UTC).isoformat(),
        "python": sys.version,
        "platform": platform.platform(),
        "cpu_model": cpu_model,
        "packages": packages,
        "ggann_source": _ggann_source_metadata(),
        "git_revision": args.git_revision,
        "ggann_backend": args.ggann_backend,
        "thread_settings": {
            name: os.environ.get(name)
            for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS")
        },
        "profiling_warning": (
            "cProfile and tracemalloc timings are diagnostic and intentionally separate "
            "from the uninstrumented samples and matched Scanpy benchmark."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preset", choices=("smoke", "standard", "extended"), default="extended")
    parser.add_argument("--format", choices=("dense", "csr", "csc"), default="csr")
    parser.add_argument("--workloads", default=",".join(_WORKLOADS))
    parser.add_argument("--stages", default=",".join(_STAGES))
    parser.add_argument("--source", choices=("x", "layer", "raw"), default="x")
    parser.add_argument(
        "--ggann-backend",
        choices=("plotnine", "matplotlib"),
        default="plotnine",
    )
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--seed", type=int, default=20260809)
    parser.add_argument("--rss-interval-ms", type=float, default=1.0)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--profile-dir", type=Path, required=True)
    parser.add_argument("--git-revision", default="unknown")
    args = parser.parse_args()

    workloads = [item.strip() for item in args.workloads.split(",") if item.strip()]
    stages = [item.strip() for item in args.stages.split(",") if item.strip()]
    invalid_workloads = sorted(set(workloads) - set(_WORKLOADS))
    invalid_stages = sorted(set(stages) - set(_STAGES))
    if invalid_workloads:
        parser.error(f"unknown workload(s): {', '.join(invalid_workloads)}")
    if invalid_stages:
        parser.error(f"unknown stage(s): {', '.join(invalid_stages)}")
    if args.repeats < 1:
        parser.error("--repeats must be positive")

    shape = _shape_for(args.preset, args.format)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.profile_dir.mkdir(parents=True, exist_ok=True)
    results = []
    with tempfile.TemporaryDirectory(prefix="ggann-profile-") as temporary:
        output_dir = Path(temporary)
        for workload in workloads:
            spec = ComparisonSpec(
                preset=args.preset,
                matrix_format=args.format,
                workload=workload,
                source=args.source,
                repeats=args.repeats,
                seed=args.seed,
                rss_interval_seconds=args.rss_interval_ms / 1000,
                ggann_backend=args.ggann_backend,
                **shape,
            )
            adata, genes, input_bytes = _make_fixture(spec.fixture_spec())
            _add_fixture_metadata(adata)
            before = _adata_digest(adata)

            # Warm plotnine, Matplotlib, and font caches before measured/profiled calls.
            warm = _construct("ggann", adata, genes, spec)
            warm_figure = warm.draw(show=False)
            from matplotlib import pyplot as plt

            plt.close(warm_figure)
            _close_output(warm)
            del warm, warm_figure
            gc.collect()

            stage_results = [
                _measure_stage(stage, adata, genes, spec, output_dir, args.profile_dir)
                for stage in stages
            ]
            after = _adata_digest(adata)
            results.append(
                {
                    "case_id": spec.case_id,
                    "parameters": asdict(spec),
                    "input_bytes": input_bytes,
                    "input_fingerprint_before": before,
                    "input_fingerprint_after": after,
                    "input_immutable": before == after,
                    "stages": stage_results,
                }
            )

    payload = {
        "schema_version": 1,
        "profile_kind": "ggann_hotspots",
        "metadata": _metadata(args),
        "results": results,
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
