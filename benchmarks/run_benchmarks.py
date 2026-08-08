#!/usr/bin/env python3
"""Run deterministic ggann preparation and rendering benchmarks.

Each case runs in a fresh child process. Fixture construction and imports happen
before measurement, so ``cold`` means the first workload call in that process.
"""

from __future__ import annotations

import argparse
import base64
import gc
import hashlib
import importlib.metadata
import json
import os
import platform
import statistics
import subprocess
import sys
import tempfile
import threading
import time
from dataclasses import asdict, dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, Callable


for _name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_name, "1")
os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault(
    "MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "ggann-benchmark-mpl")
)


_RESULT_PREFIX = "GGANN_BENCH_RESULT="

_SCALING_FIELDS = ("n_obs", "n_vars", "n_genes", "n_groups")

_PRESETS = {
    "smoke": {
        "default": {
            "n_obs": 500,
            "n_vars": 200,
            "n_genes": 4,
            "n_groups": 4,
            "density": 0.03,
            "embedding_dims": 10,
            "render_cells": 300,
        }
    },
    "standard": {
        "default": {
            "n_obs": 4_000,
            "n_vars": 2_000,
            "n_genes": 16,
            "n_groups": 16,
            "density": 0.005,
            "embedding_dims": 50,
            "render_cells": 2_000,
        }
    },
    "extended": {
        "dense": {
            "n_obs": 6_000,
            "n_vars": 3_000,
            "n_genes": 32,
            "n_groups": 24,
            "density": 0.003,
            "embedding_dims": 50,
            "render_cells": 3_000,
        },
        "sparse": {
            "n_obs": 20_000,
            "n_vars": 10_000,
            "n_genes": 32,
            "n_groups": 32,
            "density": 0.001,
            "embedding_dims": 50,
            "render_cells": 5_000,
        },
    },
}

_CORE_WORKLOADS = (
    "resolve_x",
    "tidy_x",
    "aggregate_x",
)

_PREPARATION_WORKLOADS = (
    *_CORE_WORKLOADS,
    "resolve_layer",
    "resolve_raw",
    "resolve_mixed",
    "resolve_pca",
    "aggregate_layer",
    "aggregate_raw",
    "group_means_x",
    "grammar_prepare",
    "plot_embedding_prepare",
    "plot_features_prepare",
    "plot_dotplot_prepare",
    "plot_highest_expr_prepare",
)

_RENDER_WORKLOADS = (
    "render_embedding",
    "render_dotplot",
)


@dataclass(frozen=True)
class CaseSpec:
    preset: str
    matrix_format: str
    workload: str
    n_obs: int
    n_vars: int
    n_genes: int
    n_groups: int
    density: float
    embedding_dims: int
    render_cells: int
    seed: int
    repeats: int
    rss_interval_seconds: float

    @property
    def case_id(self) -> str:
        defaults = _shape_for(self.preset, self.matrix_format)
        changed_shape = [
            f"{name}={getattr(self, name)}"
            for name in _SCALING_FIELDS
            if getattr(self, name) != defaults[name]
        ]
        shape_suffix = f"[{','.join(changed_shape)}]" if changed_shape else ""
        return f"{self.preset}{shape_suffix}/{self.matrix_format}/{self.workload}"


class _RSSSampler:
    def __init__(self, interval_seconds: float):
        self.interval_seconds = interval_seconds
        self.peak_bytes = _rss_bytes()
        self._stop = threading.Event()
        self._thread = threading.Thread(target=self._sample, daemon=True)

    def _sample(self) -> None:
        while not self._stop.wait(self.interval_seconds):
            self.peak_bytes = max(self.peak_bytes, _rss_bytes())

    def __enter__(self) -> _RSSSampler:
        self._thread.start()
        return self

    def __exit__(self, *_exc: object) -> None:
        self.peak_bytes = max(self.peak_bytes, _rss_bytes())
        self._stop.set()
        self._thread.join()


def _rss_bytes() -> int:
    """Current resident set size, using Linux procfs when available."""
    statm = Path("/proc/self/statm")
    if statm.exists():
        resident_pages = int(statm.read_text().split()[1])
        return resident_pages * os.sysconf("SC_PAGE_SIZE")

    # Fallback is a process high-water mark, not current RSS. It keeps the runner
    # usable off Linux, but retained-memory measurements then equal the high-water
    # mark and should not be compared with Linux results.
    import resource

    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(rss if sys.platform == "darwin" else rss * 1024)


def _matrix_bytes(matrix: Any) -> int:
    from scipy import sparse

    if sparse.issparse(matrix):
        return int(matrix.data.nbytes + matrix.indices.nbytes + matrix.indptr.nbytes)
    return int(getattr(matrix, "nbytes", 0))


def _frame_description(frame: Any, *, fingerprint: bool) -> dict[str, Any]:
    import pandas as pd

    frame_bytes = int(frame.memory_usage(index=True, deep=True).sum())
    categorical = {}
    for column in frame.columns:
        series = frame[column]
        if isinstance(series.dtype, pd.CategoricalDtype):
            categorical[str(column)] = {
                "categories": [str(value) for value in series.cat.categories],
                "ordered": bool(series.cat.ordered),
            }
    description = {
        "kind": "dataframe",
        "shape": [int(frame.shape[0]), int(frame.shape[1])],
        "bytes": frame_bytes,
        "columns": [str(column) for column in frame.columns],
        "dtypes": [str(dtype) for dtype in frame.dtypes],
        "categorical": categorical,
    }
    if fingerprint:
        try:
            hashed = pd.util.hash_pandas_object(frame, index=True, categorize=True)
            digest = hashlib.sha256()
            digest.update(json.dumps(description, sort_keys=True).encode())
            digest.update(hashed.to_numpy().tobytes())
            description["fingerprint"] = digest.hexdigest()
        except (TypeError, ValueError):
            description["fingerprint"] = hashlib.sha256(
                repr(frame).encode()
            ).hexdigest()
    return description


def _output_description(output: Any, *, fingerprint: bool) -> dict[str, Any]:
    import numpy as np
    import pandas as pd

    if isinstance(output, pd.DataFrame):
        return _frame_description(output, fingerprint=fingerprint)
    plot_data = getattr(output, "data", None)
    if isinstance(plot_data, pd.DataFrame):
        description = _frame_description(plot_data, fingerprint=fingerprint)
        description["kind"] = type(output).__name__
        return description
    if isinstance(output, np.ndarray):
        description = {
            "kind": "ndarray",
            "shape": [int(value) for value in output.shape],
            "bytes": int(output.nbytes),
        }
        if fingerprint:
            description["fingerprint"] = hashlib.sha256(output.tobytes()).hexdigest()
        return description

    encoded = json.dumps(output, sort_keys=True, default=str).encode()
    description = {
        "kind": type(output).__name__,
        "shape": None,
        "bytes": len(encoded),
    }
    if fingerprint:
        description["fingerprint"] = hashlib.sha256(encoded).hexdigest()
    return description


def _matrix_description(matrix: Any) -> dict[str, Any]:
    """Describe matrix storage without converting its representation."""
    from scipy import sparse

    description = {
        "kind": f"sparse_{matrix.format}" if sparse.issparse(matrix) else "ndarray",
        "shape": [int(value) for value in matrix.shape],
        "bytes": _matrix_bytes(matrix),
        "dtype": str(matrix.dtype),
    }
    if sparse.issparse(matrix):
        description["nnz"] = int(matrix.nnz)
    return description


def _stage_sizes(
    adata: Any,
    genes: list[str],
    workload: str,
    output_description: dict[str, Any],
) -> dict[str, dict[str, Any]]:
    """Measure major preparation stages after all timed workload calls."""
    prepared_fields = ("kind", "shape", "bytes", "columns", "dtypes")
    stages = {
        "prepared_data": {
            name: output_description[name]
            for name in prepared_fields
            if name in output_description
        }
    }
    if workload in _CORE_WORKLOADS:
        stages["projected_expression"] = _matrix_description(adata[:, genes].X)
        stages["projected_obs"] = _frame_description(
            adata.obs.loc[:, ["group", "split"]], fingerprint=False
        )
    elif workload == "plot_highest_expr_prepare":
        stages["whole_expression"] = _matrix_description(adata.X)
    return stages


def _measure_call(
    function: Callable[[], Any], rss_interval_seconds: float
) -> dict[str, Any]:
    gc.collect()
    baseline_rss = _rss_bytes()
    with _RSSSampler(rss_interval_seconds) as sampler:
        start = time.perf_counter_ns()
        output = function()
        duration_ns = time.perf_counter_ns() - start

    output_rss = _rss_bytes()
    peak_rss = max(sampler.peak_bytes, output_rss)
    description = _output_description(output, fingerprint=False)
    del output
    gc.collect()
    final_rss = _rss_bytes()
    return {
        "duration_seconds": duration_ns / 1_000_000_000,
        "baseline_rss_bytes": baseline_rss,
        "peak_rss_bytes": peak_rss,
        "peak_rss_delta_bytes": max(0, peak_rss - baseline_rss),
        "retained_with_output_bytes": max(0, output_rss - baseline_rss),
        "retained_after_gc_bytes": max(0, final_rss - baseline_rss),
        "output": description,
    }


def _copy_matrix(matrix: Any) -> Any:
    return matrix.copy()


def _transform_matrix(matrix: Any, function: Callable[[Any], Any]) -> Any:
    import numpy as np
    from scipy import sparse

    if sparse.issparse(matrix):
        result = matrix.copy()
        result.data = np.asarray(function(result.data), dtype=np.float32)
        return result
    return np.asarray(function(matrix), dtype=np.float32)


def _make_fixture(spec: CaseSpec) -> tuple[Any, list[str], dict[str, int]]:
    import numpy as np
    import pandas as pd
    from anndata import AnnData
    from scipy import sparse

    rng = np.random.default_rng(spec.seed)
    counts_sparse = sparse.random(
        spec.n_obs,
        spec.n_vars,
        density=spec.density,
        format="csr",
        dtype=np.float32,
        random_state=rng,
        data_rvs=lambda size: rng.integers(1, 8, size=size).astype(np.float32),
    )
    if spec.matrix_format == "dense":
        counts = counts_sparse.toarray()
    elif spec.matrix_format == "csr":
        counts = counts_sparse
    elif spec.matrix_format == "csc":
        counts = counts_sparse.tocsc()
    else:  # guarded by argument validation
        raise ValueError(f"Unknown matrix format: {spec.matrix_format}")

    obs_names = pd.Index([f"cell_{index}" for index in range(spec.n_obs)])
    var_names = pd.Index([f"gene_{index}" for index in range(spec.n_vars)])
    group_names = [f"group_{index}" for index in range(spec.n_groups)]
    group_codes = np.arange(spec.n_obs) % spec.n_groups
    rng.shuffle(group_codes)
    split_codes = np.arange(spec.n_obs) % 2
    rng.shuffle(split_codes)
    obs = pd.DataFrame(
        {
            "group": pd.Categorical.from_codes(group_codes, categories=group_names),
            "split": pd.Categorical.from_codes(split_codes, categories=["A", "B"]),
            "score": rng.normal(size=spec.n_obs).astype(np.float32),
        },
        index=obs_names,
    )
    var = pd.DataFrame(index=var_names)

    x_matrix = _transform_matrix(counts, np.log1p)
    adata = AnnData(X=x_matrix, obs=obs, var=var)
    adata.layers["counts"] = _copy_matrix(counts)
    adata.layers["sqrt_counts"] = _transform_matrix(counts, np.sqrt)
    raw = AnnData(
        X=_copy_matrix(counts),
        obs=pd.DataFrame(index=obs_names.copy()),
        var=pd.DataFrame(index=var_names.copy()),
    )
    adata.raw = raw
    adata.obsm["X_umap"] = rng.normal(size=(spec.n_obs, 2)).astype(np.float32)
    adata.obsm["X_pca"] = rng.normal(size=(spec.n_obs, spec.embedding_dims)).astype(
        np.float32
    )

    selected_indices = np.linspace(
        0, spec.n_vars - 1, num=min(spec.n_genes, spec.n_vars), dtype=int
    )
    selected_genes = [
        str(var_names[index]) for index in dict.fromkeys(selected_indices)
    ]

    input_bytes = {
        "x": _matrix_bytes(adata.X),
        "layers": sum(_matrix_bytes(matrix) for matrix in adata.layers.values()),
        "raw": _matrix_bytes(adata.raw.X),
        "obsm": sum(_matrix_bytes(matrix) for matrix in adata.obsm.values()),
        "obs": int(adata.obs.memory_usage(index=True, deep=True).sum()),
        "var": int(adata.var.memory_usage(index=True, deep=True).sum()),
    }
    input_bytes["total"] = sum(input_bytes.values())
    return adata, selected_genes, input_bytes


def _workload_function(
    adata: Any, genes: list[str], name: str, render_cells: int
) -> Callable[[], Any]:
    if name == "resolve_x":
        from ggann._resolve import resolve_frame

        return lambda: resolve_frame(adata, ["group", "split", *genes], use_raw=False)
    if name == "resolve_layer":
        from ggann._resolve import resolve_frame

        return lambda: resolve_frame(adata, ["group", "split", *genes], layer="counts")
    if name == "resolve_raw":
        from ggann._resolve import resolve_frame

        return lambda: resolve_frame(adata, ["group", "split", *genes], use_raw=True)
    if name == "resolve_mixed":
        from ggann import gene, obs, obsm
        from ggann._resolve import resolve_frame

        refs = [
            obs("group"),
            gene(genes[0], use_raw=False),
            gene(genes[1], layer="counts"),
            gene(genes[2], use_raw=True),
            obsm("umap", 0),
            obsm("umap", 1),
        ]
        return lambda: resolve_frame(adata, refs, use_raw=False)
    if name == "resolve_pca":
        from ggann._resolve import embedding_coords

        return lambda: embedding_coords(adata, "pca", n=2)
    if name == "tidy_x":
        from ggann._aggregate import tidy_expression

        return lambda: tidy_expression(
            adata,
            genes,
            "group",
            use_raw=False,
            extra_obs=("split",),
        )
    if name == "aggregate_x":
        from ggann._aggregate import aggregate_expression

        return lambda: aggregate_expression(
            adata,
            genes,
            "group",
            use_raw=False,
            extra_by=("split",),
        )
    if name == "aggregate_layer":
        from ggann._aggregate import aggregate_expression

        return lambda: aggregate_expression(
            adata,
            genes,
            "group",
            layer="counts",
            extra_by=("split",),
        )
    if name == "aggregate_raw":
        from ggann._aggregate import aggregate_expression

        return lambda: aggregate_expression(
            adata,
            genes,
            "group",
            use_raw=True,
            extra_by=("split",),
        )
    if name == "group_means_x":
        from ggann._aggregate import group_means

        return lambda: group_means(adata, genes, "group", use_raw=False)
    if name == "grammar_prepare":
        from ggann import aes, gganndata

        return lambda: gganndata(
            adata,
            aes("UMAP_1", "UMAP_2", color=genes[0]),
            use_raw=False,
        )
    if name == "plot_embedding_prepare":
        from ggann import plot_embedding

        return lambda: plot_embedding(
            adata,
            color="group",
            pointdensity=False,
            downsample=min(render_cells, adata.n_obs),
        )
    if name == "plot_features_prepare":
        from ggann import plot_features

        return lambda: plot_features(
            adata,
            genes[: min(4, len(genes))],
            use_raw=False,
            downsample=min(render_cells, adata.n_obs),
        )
    if name == "plot_dotplot_prepare":
        from ggann import plot_dotplot

        return lambda: plot_dotplot(
            adata, genes, "group", use_raw=False, split_by="split"
        )
    if name == "plot_highest_expr_prepare":
        from ggann import plot_highest_expr_genes

        return lambda: plot_highest_expr_genes(
            adata, n=min(20, adata.n_vars), use_raw=False
        )
    if name == "render_embedding":
        from ggann import plot_embedding

        def render_embedding() -> dict[str, int]:
            from matplotlib import pyplot as plt

            plot = plot_embedding(
                adata,
                color="group",
                pointdensity=False,
                downsample=min(render_cells, adata.n_obs),
            )
            figure = plot.draw(show=False)
            summary = {"axes": len(figure.axes), "rows": len(plot.data)}
            plt.close(figure)
            return summary

        return render_embedding
    if name == "render_dotplot":
        from ggann import plot_dotplot

        def render_dotplot() -> dict[str, int]:
            from matplotlib import pyplot as plt

            plot = plot_dotplot(adata, genes, "group", use_raw=False, split_by="split")
            figure = plot.draw(show=False)
            summary = {"axes": len(figure.axes), "rows": len(plot.data)}
            plt.close(figure)
            return summary

        return render_dotplot
    raise ValueError(f"Unknown workload: {name}")


def _execute_case(spec: CaseSpec) -> dict[str, Any]:
    adata, genes, input_bytes = _make_fixture(spec)
    function = _workload_function(adata, genes, spec.workload, spec.render_cells)

    cold = _measure_call(function, spec.rss_interval_seconds)
    repeated_samples = [
        _measure_call(function, spec.rss_interval_seconds) for _ in range(spec.repeats)
    ]
    durations = [sample["duration_seconds"] for sample in repeated_samples]
    measured_outputs = [
        cold["output"],
        *(sample["output"] for sample in repeated_samples),
    ]
    semantic_output = function()
    output = _output_description(semantic_output, fingerprint=True)
    stage_sizes = _stage_sizes(adata, genes, spec.workload, output)
    del semantic_output
    gc.collect()
    repeated = {
        "samples": repeated_samples,
        "median_duration_seconds": statistics.median(durations),
        "min_duration_seconds": min(durations),
        "max_duration_seconds": max(durations),
        "max_peak_rss_delta_bytes": max(
            sample["peak_rss_delta_bytes"] for sample in repeated_samples
        ),
        "max_retained_after_gc_bytes": max(
            sample["retained_after_gc_bytes"] for sample in repeated_samples
        ),
    }
    return {
        "case_id": spec.case_id,
        "parameters": asdict(spec),
        "selected_genes": genes,
        "input_bytes": input_bytes,
        "cold": cold,
        "repeated": repeated,
        "output": output,
        "stage_sizes": stage_sizes,
        "output_stable": all(item == measured_outputs[0] for item in measured_outputs),
    }


def _encoded_case(spec: CaseSpec) -> str:
    encoded = json.dumps(asdict(spec), separators=(",", ":")).encode()
    return base64.urlsafe_b64encode(encoded).decode()


def _decoded_case(payload: str) -> CaseSpec:
    decoded = base64.urlsafe_b64decode(payload.encode())
    return CaseSpec(**json.loads(decoded))


def _run_case_process(spec: CaseSpec, timeout_seconds: float) -> dict[str, Any]:
    environment = os.environ.copy()
    environment["PYTHONHASHSEED"] = "0"
    environment["MPLBACKEND"] = "Agg"
    command = [
        sys.executable,
        str(Path(__file__).resolve()),
        "--internal-case",
        _encoded_case(spec),
    ]
    completed = subprocess.run(
        command,
        check=False,
        capture_output=True,
        text=True,
        env=environment,
        timeout=timeout_seconds,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            f"Benchmark case {spec.case_id} failed with exit code "
            f"{completed.returncode}:\n{completed.stderr.strip()}"
        )
    result_lines = [
        line.removeprefix(_RESULT_PREFIX)
        for line in completed.stdout.splitlines()
        if line.startswith(_RESULT_PREFIX)
    ]
    if len(result_lines) != 1:
        raise RuntimeError(
            f"Benchmark case {spec.case_id} returned no parseable result.\n"
            f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
        )
    result = json.loads(result_lines[0])
    diagnostics = [line for line in completed.stderr.splitlines() if line.strip()]
    if diagnostics:
        result["diagnostics"] = diagnostics
    return result


def _shape_for(preset: str, matrix_format: str) -> dict[str, Any]:
    config = _PRESETS[preset]
    if preset != "extended":
        return dict(config["default"])
    key = "dense" if matrix_format == "dense" else "sparse"
    return dict(config[key])


def _positive_int_values(value: str) -> list[int]:
    """Parse one or more comma-separated positive integers."""
    try:
        values = [int(item.strip()) for item in value.split(",") if item.strip()]
    except ValueError as error:
        raise argparse.ArgumentTypeError(
            "expected a positive integer or comma-separated positive integers"
        ) from error
    if not values or any(item < 1 for item in values):
        raise argparse.ArgumentTypeError(
            "expected a positive integer or comma-separated positive integers"
        )
    return list(dict.fromkeys(values))


def _shape_variants(
    preset: str,
    matrix_format: str,
    overrides: dict[str, list[int] | None],
) -> list[dict[str, Any]]:
    """Return shapes varying at most one dimension within a benchmark run."""
    provided = {
        name: values for name, values in overrides.items() if values is not None
    }
    varying = [name for name, values in provided.items() if len(values) > 1]
    if len(varying) > 1:
        raise ValueError(
            "Only one of --n-obs, --n-vars, --n-genes and --n-groups may "
            "contain multiple values in one run."
        )

    base = _shape_for(preset, matrix_format)
    scale_field = varying[0] if varying else None
    scale_values = provided.get(scale_field, [None])
    shapes = []
    for scale_value in scale_values:
        shape = dict(base)
        for name, values in provided.items():
            shape[name] = scale_value if name == scale_field else values[0]
        if shape["n_genes"] > shape["n_vars"]:
            raise ValueError(
                f"n_genes ({shape['n_genes']}) cannot exceed n_vars "
                f"({shape['n_vars']})."
            )
        if shape["n_groups"] > shape["n_obs"]:
            raise ValueError(
                f"n_groups ({shape['n_groups']}) cannot exceed n_obs "
                f"({shape['n_obs']})."
            )
        shapes.append(shape)
    return shapes


def _selected_formats(value: str) -> list[str]:
    formats = [item.strip().lower() for item in value.split(",") if item.strip()]
    if not formats:
        raise ValueError("At least one matrix format is required.")
    if formats == ["all"]:
        formats = ["dense", "csr", "csc"]
    invalid = sorted(set(formats) - {"dense", "csr", "csc"})
    if invalid:
        raise ValueError(f"Unknown matrix format(s): {', '.join(invalid)}")
    return list(dict.fromkeys(formats))


def _selected_workloads(value: str, include_rendering: bool) -> list[str]:
    if value == "core":
        workloads = list(_CORE_WORKLOADS)
    elif value == "all":
        workloads = list(_PREPARATION_WORKLOADS)
    else:
        workloads = [item.strip() for item in value.split(",") if item.strip()]
    if include_rendering:
        workloads.extend(_RENDER_WORKLOADS)
    if not workloads:
        raise ValueError("At least one workload is required.")
    available = set(_PREPARATION_WORKLOADS) | set(_RENDER_WORKLOADS)
    invalid = sorted(set(workloads) - available)
    if invalid:
        raise ValueError(f"Unknown workload(s): {', '.join(invalid)}")
    return list(dict.fromkeys(workloads))


def _git_metadata() -> dict[str, Any]:
    def git(*arguments: str) -> str | None:
        completed = subprocess.run(
            ["git", *arguments], check=False, capture_output=True, text=True
        )
        return completed.stdout.strip() if completed.returncode == 0 else None

    status = git("status", "--porcelain")
    return {
        "revision": git("rev-parse", "HEAD"),
        "dirty": bool(status) if status is not None else None,
    }


def _versions() -> dict[str, str | None]:
    versions = {}
    for package in (
        "ggann",
        "annplyr",
        "anndata",
        "numpy",
        "pandas",
        "scipy",
        "plotnine",
        "plotnine-extra",
    ):
        try:
            versions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            versions[package] = None
    return versions


def _ggann_source_metadata() -> dict[str, Any]:
    """Describe the imported ggann tree used by benchmark child processes."""
    import ggann

    package_file = Path(ggann.__file__).resolve()
    package_root = package_file.parent
    source_files = sorted(path for path in package_root.rglob("*.py") if path.is_file())
    digest = hashlib.sha256()
    for path in source_files:
        relative = path.relative_to(package_root).as_posix().encode()
        digest.update(len(relative).to_bytes(8, "big"))
        digest.update(relative)
        contents = path.read_bytes()
        digest.update(len(contents).to_bytes(8, "big"))
        digest.update(contents)
    return {
        "package_file": str(package_file),
        "package_root": str(package_root),
        "python_files": len(source_files),
        "python_tree_sha256": digest.hexdigest(),
    }


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preset", choices=tuple(_PRESETS), default="standard")
    parser.add_argument(
        "--formats",
        default="dense,csr,csc",
        help="Comma-separated dense, csr, csc, or all.",
    )
    parser.add_argument(
        "--workloads",
        default="core",
        help="core, all, or a comma-separated workload list.",
    )
    parser.add_argument(
        "--include-rendering",
        action="store_true",
        help="Add bounded embedding and dotplot rendering cases.",
    )
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--seed", type=int, default=20_260_808)
    parser.add_argument("--rss-interval-ms", type=float, default=1.0)
    parser.add_argument("--timeout-seconds", type=float, default=900.0)
    parser.add_argument("--label", default="benchmark")
    parser.add_argument("--output", type=Path)
    for field in _SCALING_FIELDS:
        option = "--" + field.replace("_", "-")
        parser.add_argument(
            option,
            type=_positive_int_values,
            help=(
                "Override the preset with one value, or provide comma-separated "
                "values to scale this dimension. Only one dimension may vary per run."
            ),
        )
    parser.add_argument(
        "--list-workloads", action="store_true", help="List workload names and exit."
    )
    parser.add_argument("--internal-case", help=argparse.SUPPRESS)
    return parser


def _internal_main(payload: str) -> int:
    try:
        result = _execute_case(_decoded_case(payload))
    except Exception:
        import traceback

        traceback.print_exc()
        return 1
    print(_RESULT_PREFIX + json.dumps(result, separators=(",", ":")))
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    if args.internal_case:
        return _internal_main(args.internal_case)
    if args.list_workloads:
        for workload in (*_PREPARATION_WORKLOADS, *_RENDER_WORKLOADS):
            print(workload)
        return 0
    if args.repeats < 1:
        parser.error("--repeats must be at least 1")
    if args.rss_interval_ms <= 0:
        parser.error("--rss-interval-ms must be positive")
    if args.output is None:
        parser.error("--output is required")

    try:
        formats = _selected_formats(args.formats)
        workloads = _selected_workloads(args.workloads, args.include_rendering)
        shape_overrides = {field: getattr(args, field) for field in _SCALING_FIELDS}
        shapes_by_format = {
            matrix_format: _shape_variants(args.preset, matrix_format, shape_overrides)
            for matrix_format in formats
        }
    except ValueError as error:
        parser.error(str(error))

    cases = []
    for matrix_format in formats:
        for shape in shapes_by_format[matrix_format]:
            for workload in workloads:
                cases.append(
                    CaseSpec(
                        preset=args.preset,
                        matrix_format=matrix_format,
                        workload=workload,
                        seed=args.seed,
                        repeats=args.repeats,
                        rss_interval_seconds=args.rss_interval_ms / 1_000,
                        **shape,
                    )
                )

    # Capture provenance before any child starts so the digest describes the
    # source tree used for the run even if a long benchmark overlaps later edits.
    package_versions = _versions()
    ggann_source = _ggann_source_metadata()

    results = []
    for index, spec in enumerate(cases, start=1):
        print(f"[{index}/{len(cases)}] {spec.case_id}", flush=True)
        results.append(_run_case_process(spec, args.timeout_seconds))

    document = {
        "schema_version": 1,
        "metadata": {
            "label": args.label,
            "created_at": datetime.now(UTC).isoformat(),
            "python": sys.version,
            "platform": platform.platform(),
            "rss_backend": "proc_statm"
            if Path("/proc/self/statm").exists()
            else "resource_peak",
            "packages": package_versions,
            "ggann_source": ggann_source,
            "thread_settings": {
                name: os.environ.get(name)
                for name in (
                    "OPENBLAS_NUM_THREADS",
                    "OMP_NUM_THREADS",
                    "MKL_NUM_THREADS",
                )
            },
            "git": _git_metadata(),
            "preset": args.preset,
            "formats": formats,
            "workloads": workloads,
            "shape_overrides": {
                name: values
                for name, values in shape_overrides.items()
                if values is not None
            },
            "repeats": args.repeats,
            "seed": args.seed,
            "rss_interval_ms": args.rss_interval_ms,
        },
        "results": results,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(document, indent=2) + "\n")
    print(f"Wrote {len(results)} cases to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
