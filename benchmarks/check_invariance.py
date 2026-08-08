#!/usr/bin/env python3
"""Snapshot and compare representative ggann downstream artifacts."""

from __future__ import annotations

import argparse
import importlib.metadata
import json
import math
import os
from pathlib import Path
import platform
import sys
from typing import Any

import numpy as np
import pandas as pd


def _mapping(plot) -> dict[str, str]:
    return {str(key): str(value) for key, value in plot.mapping.items()}


def _json_value(value):
    """Return a stable JSON representation for public plot scale properties."""
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, (float, np.floating)):
        return value if math.isfinite(value) else str(value)
    if isinstance(value, (list, tuple)):
        return [_json_value(item) for item in value]
    if callable(value):
        samples = [0.0, 0.25, 0.5, 0.75, 1.0]
        try:
            return {"callable_samples": _json_value(value(samples))}
        except (TypeError, ValueError):
            return {
                "callable": getattr(value, "__qualname__", type(value).__qualname__)
            }
    return str(value)


def _plot_contract(plot) -> dict[str, Any]:
    """Describe public mappings, labels and scales without object addresses."""
    scales = []
    for scale in plot.scales:
        spec = {
            "class": f"{type(scale).__module__}.{type(scale).__name__}",
            "aesthetics": list(scale.aesthetics),
        }
        for attribute in (
            "name",
            "breaks",
            "limits",
            "labels",
            "guide",
            "na_value",
            "trans",
            "range",
        ):
            try:
                spec[attribute] = _json_value(getattr(scale, attribute))
            except AttributeError:
                continue
        scales.append(spec)
    return {
        "mapping": _mapping(plot),
        "labels": {key: _json_value(value) for key, value in vars(plot.labels).items()},
        "scales": scales,
    }


def _provenance() -> dict[str, Any]:
    try:
        from benchmarks.run_benchmarks import _ggann_source_metadata
    except ModuleNotFoundError:  # direct ``python benchmarks/check_invariance.py``
        from run_benchmarks import _ggann_source_metadata

    packages = {}
    for package in (
        "ggann",
        "annplyr",
        "anndata",
        "numpy",
        "pandas",
        "scipy",
        "plotnine",
        "plotnine-extra",
        "scanpy",
        "matplotlib",
        "pillow",
    ):
        try:
            packages[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            packages[package] = None
    return {
        "python": sys.version,
        "platform": platform.platform(),
        "packages": packages,
        "thread_settings": {
            name: os.environ.get(name)
            for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS")
        },
        "ggann_source": _ggann_source_metadata(),
    }


def snapshot(output: Path) -> None:
    import scanpy as sc

    import ggann as ag
    from ggann._aggregate import aggregate_expression
    from ggann._resolve import resolve_frame

    output.mkdir(parents=True, exist_ok=True)
    adata = sc.datasets.pbmc68k_reduced()
    genes = [
        gene for gene in ("CD3D", "NKG7", "CST3", "GNLY") if gene in adata.raw.var_names
    ]
    group_by = "bulk_labels"

    plots = {
        "embedding": ag.plot_embedding(
            adata, "umap", color=genes[0], use_raw=True, downsample=300
        ),
        "dotplot": ag.plot_dotplot(adata, genes, group_by),
        "matrixplot": ag.plot_matrixplot(adata, genes, group_by),
        "heatmap": ag.plot_heatmap(
            adata,
            genes,
            group_by,
            downsample=40,
        ),
    }
    frames = {
        "resolve": resolve_frame(
            adata,
            [group_by, genes[0], ag.obsm("umap", 0), ag.obsm("umap", 1)],
        ),
        "aggregate": aggregate_expression(adata, genes, group_by),
        **{f"plot_{name}": plot.data for name, plot in plots.items()},
    }

    manifest: dict[str, Any] = {
        "schema_version": 2,
        "provenance": _provenance(),
        "dataset": "scanpy.datasets.pbmc68k_reduced",
        "shape": list(adata.shape),
        "genes": genes,
        "group_by": group_by,
        "parameters": {
            "embedding_downsample": 300,
            "heatmap_downsample_per_group": 40,
            "random_state": 0,
            "render_width": 6,
            "render_height": 4.5,
            "render_dpi": 80,
        },
        "frames": {},
        "mappings": {name: _mapping(plot) for name, plot in plots.items()},
        "plot_contracts": {name: _plot_contract(plot) for name, plot in plots.items()},
        "images": {},
    }
    for name, frame in frames.items():
        filename = f"{name}.pkl"
        frame.to_pickle(output / filename)
        manifest["frames"][name] = {
            "file": filename,
            "shape": list(frame.shape),
            "columns": [str(column) for column in frame.columns],
        }

    for name, plot in plots.items():
        filename = f"{name}.png"
        plot.save(output / filename, width=6, height=4.5, dpi=80, verbose=False)
        manifest["images"][name] = filename

    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


def _image(path: Path) -> np.ndarray:
    from PIL import Image

    return np.asarray(Image.open(path).convert("RGBA"), dtype=np.float32) / 255.0


def compare(baseline: Path, candidate: Path, report: Path | None) -> bool:
    before = json.loads((baseline / "manifest.json").read_text())
    after = json.loads((candidate / "manifest.json").read_text())
    failures: list[str] = []
    lines = [
        "# Downstream artifact invariance",
        "",
        "| Artifact | Result | Detail |",
        "|---|:---:|---|",
    ]

    comparable = True
    for name in (
        "schema_version",
        "dataset",
        "shape",
        "genes",
        "group_by",
        "parameters",
    ):
        if before.get(name) != after.get(name):
            failures.append(f"manifest {name} differs")
            comparable = False
    before_provenance = before.get("provenance", {})
    after_provenance = after.get("provenance", {})
    for name in ("python", "platform", "packages", "thread_settings"):
        if before_provenance.get(name) != after_provenance.get(name):
            failures.append(f"provenance {name} differs")
            comparable = False
    lines.append(
        f"| provenance and parameters | {'pass' if comparable else 'fail'} | "
        "exact comparison; source digests recorded separately |"
    )

    for name, spec in before["frames"].items():
        left = pd.read_pickle(baseline / spec["file"])
        right = pd.read_pickle(candidate / after["frames"][name]["file"])
        detail = "equal within rtol 1e-6, atol 1e-7"
        if name == "plot_matrixplot":
            removable = [
                column for column in left.columns if column not in right.columns
            ]
            if removable == ["fraction"]:
                left = left[right.columns]
                detail = "equal; unused baseline fraction column removed"
        try:
            pd.testing.assert_frame_equal(
                left,
                right,
                check_dtype=False,
                check_exact=False,
                check_categorical=True,
                rtol=1e-6,
                atol=1e-7,
            )
        except AssertionError as exc:
            failures.append(f"frame {name}: {exc}")
            lines.append(f"| `{name}` data | fail | values or structure changed |")
        else:
            lines.append(f"| `{name}` data | pass | {detail} |")

    mappings_equal = before["mappings"] == after["mappings"]
    if not mappings_equal:
        failures.append("plot mappings changed")
    lines.append(
        f"| plot mappings | {'pass' if mappings_equal else 'fail'} | exact comparison |"
    )

    contracts_equal = before.get("plot_contracts") == after.get("plot_contracts")
    if not contracts_equal:
        failures.append("plot labels or scales changed")
    lines.append(
        f"| plot labels and scales | {'pass' if contracts_equal else 'fail'} | "
        "exact comparison |"
    )

    for name, filename in before["images"].items():
        left = _image(baseline / filename)
        right = _image(candidate / after["images"][name])
        if left.shape != right.shape:
            failures.append(f"image {name}: shape {left.shape} != {right.shape}")
            lines.append(f"| `{name}` image | fail | image dimensions changed |")
            continue
        difference = np.abs(left - right)
        maximum = float(difference.max(initial=0.0))
        mean = float(difference.mean())
        equal = maximum <= (2 / 255) and mean <= 1e-4
        if not equal:
            failures.append(f"image {name}: max={maximum:.6g}, mean={mean:.6g}")
        lines.append(
            f"| `{name}` image | {'pass' if equal else 'fail'} | "
            f"max {maximum:.3g}; mean {mean:.3g} |"
        )

    text = "\n".join(lines) + "\n"
    print(text, end="")
    if report is not None:
        report.parent.mkdir(parents=True, exist_ok=True)
        report.write_text(text)
    if failures:
        print("\n".join(failures))
    return not failures


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    snapshot_parser = subparsers.add_parser("snapshot")
    snapshot_parser.add_argument("--output", type=Path, required=True)
    compare_parser = subparsers.add_parser("compare")
    compare_parser.add_argument("baseline", type=Path)
    compare_parser.add_argument("candidate", type=Path)
    compare_parser.add_argument("--report", type=Path)
    return parser


def main() -> int:
    args = _parser().parse_args()
    if args.command == "snapshot":
        snapshot(args.output)
        return 0
    return 0 if compare(args.baseline, args.candidate, args.report) else 1


if __name__ == "__main__":
    raise SystemExit(main())
