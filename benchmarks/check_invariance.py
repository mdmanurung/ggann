#!/usr/bin/env python3
"""Snapshot and compare representative ggann downstream artifacts."""

from __future__ import annotations

import argparse
import importlib.metadata
import json
import math
import os
import platform
import sys
from pathlib import Path
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
            return {"callable": getattr(value, "__qualname__", type(value).__qualname__)}
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


def _artist_contract(figure) -> dict[str, Any]:
    from matplotlib.collections import PathCollection, PolyCollection, QuadMesh
    from matplotlib.text import Text

    plot_axes = [axis for axis in figure.axes if not hasattr(axis, "_colorbar")]
    collections = [artist for axis in plot_axes for artist in axis.collections]
    paths = [artist for artist in collections if isinstance(artist, PathCollection)]
    meshes = [artist for artist in collections if isinstance(artist, QuadMesh)]
    polygons = [
        artist
        for artist in collections
        if isinstance(artist, PolyCollection) and not isinstance(artist, QuadMesh)
    ]
    return {
        "plot_axes": len(plot_axes),
        "guide_axes": len(figure.axes) - len(plot_axes),
        "collections": len(collections),
        "collection_classes": [type(artist).__name__ for artist in collections],
        "path_collections": len(paths),
        "quad_meshes": len(meshes),
        "represented_points": sum(len(artist.get_offsets()) for artist in paths),
        "represented_tiles": sum(int(artist.get_array().size) for artist in meshes)
        + sum(len(artist.get_paths()) for artist in polygons),
        "rasterized_collections": sum(bool(artist.get_rasterized()) for artist in collections),
        "visible_text": sorted(
            {
                artist.get_text()
                for artist in figure.findobj(match=Text)
                if artist.get_visible() and artist.get_text()
            }
        ),
    }


def snapshot(output: Path, *, backend: str = "plotnine") -> None:
    import scanpy as sc

    import ggann as ag
    from ggann._aggregate import aggregate_expression
    from ggann._resolve import resolve_frame

    try:
        from benchmarks.compare_scanpy import _adata_digest
    except ModuleNotFoundError:
        from compare_scanpy import _adata_digest

    output.mkdir(parents=True, exist_ok=True)
    adata = sc.datasets.pbmc68k_reduced()
    input_fingerprint_before = _adata_digest(adata)
    genes = [gene for gene in ("CD3D", "NKG7", "CST3", "GNLY") if gene in adata.raw.var_names]
    group_by = "bulk_labels"
    backend_kwargs = {} if backend == "plotnine" else {"backend": backend}

    plots = {
        "embedding": ag.plot_embedding(
            adata,
            "umap",
            color=genes[0],
            use_raw=True,
            downsample=300,
            **backend_kwargs,
        ),
        "dotplot": ag.plot_dotplot(adata, genes, group_by, **backend_kwargs),
        "matrixplot": ag.plot_matrixplot(adata, genes, group_by, **backend_kwargs),
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
        "schema_version": 3,
        "provenance": _provenance(),
        "dataset": "scanpy.datasets.pbmc68k_reduced",
        "shape": list(adata.shape),
        "input_fingerprint_before": input_fingerprint_before,
        "genes": genes,
        "group_by": group_by,
        "parameters": {
            "embedding_downsample": 300,
            "heatmap_downsample_per_group": 40,
            "random_state": 0,
            "render_width": 6,
            "render_height": 4.5,
            "render_dpi": 80,
            "backend": backend,
        },
        "frames": {},
        "mappings": {name: _mapping(plot) for name, plot in plots.items()},
        "plot_contracts": {name: _plot_contract(plot) for name, plot in plots.items()},
        "images": {},
        "artists": {},
        "vectors": {},
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
        figure = plot.draw(show=False)
        figure.set_size_inches(6, 4.5, forward=True)
        figure.savefig(output / filename, format="png", dpi=80)
        manifest["images"][name] = filename
        manifest["artists"][name] = _artist_contract(figure)
        if name in {"embedding", "dotplot", "matrixplot"}:
            vector_name = f"{name}.svg"
            figure.savefig(output / vector_name, format="svg")
            vector_text = (output / vector_name).read_text()
            manifest["vectors"][name] = {
                "file": vector_name,
                "contains_embedded_image": "<image" in vector_text,
                "bytes": (output / vector_name).stat().st_size,
            }

    input_fingerprint_after = _adata_digest(adata)
    manifest["input_fingerprint_after"] = input_fingerprint_after
    manifest["input_immutable"] = input_fingerprint_before == input_fingerprint_after

    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


def _image(path: Path) -> np.ndarray:
    from PIL import Image

    return np.asarray(Image.open(path).convert("RGBA"), dtype=np.float32) / 255.0


def compare(
    baseline: Path,
    candidate: Path,
    report: Path | None,
    *,
    allowed_image_differences: set[str] | None = None,
    image_mean_tolerance: float = 1e-4,
    image_psnr_min: float = math.inf,
    allow_backend_difference: bool = False,
) -> bool:
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
        "input_fingerprint_before",
        "input_fingerprint_after",
        "input_immutable",
    ):
        if before.get(name) != after.get(name):
            failures.append(f"manifest {name} differs")
            comparable = False
    if not before.get("input_immutable") or not after.get("input_immutable"):
        failures.append("AnnData fingerprint changed during snapshot")
        comparable = False
    before_parameters = dict(before.get("parameters", {}))
    after_parameters = dict(after.get("parameters", {}))
    if allow_backend_difference:
        before_parameters.pop("backend", None)
        after_parameters.pop("backend", None)
    if before_parameters != after_parameters:
        failures.append("manifest parameters differs")
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
            removable = [column for column in left.columns if column not in right.columns]
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
    lines.append(f"| plot mappings | {'pass' if mappings_equal else 'fail'} | exact comparison |")

    contracts_equal = before.get("plot_contracts") == after.get("plot_contracts")
    if not contracts_equal:
        failures.append("plot labels or scales changed")
    lines.append(
        f"| plot labels and scales | {'pass' if contracts_equal else 'fail'} | exact comparison |"
    )

    for name in before["images"]:
        before_artists = before.get("artists", {}).get(name)
        after_artists = after.get("artists", {}).get(name)
        if before_artists is None or after_artists is None:
            failures.append(f"artists {name}: contract missing")
            lines.append(f"| `{name}` artists | fail | artist contract missing |")
            continue
        semantic_fields = (
            "collections",
            "represented_points",
            "represented_tiles",
            "rasterized_collections",
            "visible_text",
        )
        fields = semantic_fields if allow_backend_difference else tuple(before_artists)
        equal = all(before_artists.get(field) == after_artists.get(field) for field in fields)
        if not equal:
            failures.append(f"artists {name}: semantic contract changed")
        lines.append(
            f"| `{name}` artists | {'pass' if equal else 'fail'} | "
            f"{after_artists['represented_points']} points; "
            f"{after_artists['represented_tiles']} tiles; "
            f"{after_artists['rasterized_collections']} rasterized collections; "
            "guide and axis text exact |"
        )

    for name in ("embedding", "dotplot", "matrixplot"):
        before_vector = before.get("vectors", {}).get(name)
        after_vector = after.get("vectors", {}).get(name)
        equal = (
            before_vector is not None
            and after_vector is not None
            and not before_vector.get("contains_embedded_image", True)
            and not after_vector.get("contains_embedded_image", True)
            and (candidate / after_vector["file"]).is_file()
        )
        if not equal:
            failures.append(f"vector {name}: embedded raster or artifact missing")
        lines.append(
            f"| `{name}` SVG | {'pass' if equal else 'fail'} | "
            "vector collections; no embedded image |"
        )

    allowed_image_differences = allowed_image_differences or set()
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
        rmse = float(np.sqrt(np.mean(np.square(left - right))))
        psnr = math.inf if rmse == 0 else 20 * math.log10(1 / rmse)
        exact_enough = maximum <= (2 / 255) and mean <= 1e-4
        explicitly_tolerated = (
            name in allowed_image_differences
            and mean <= image_mean_tolerance
            and psnr >= image_psnr_min
        )
        equal = exact_enough or explicitly_tolerated
        if not equal:
            failures.append(f"image {name}: max={maximum:.6g}, mean={mean:.6g}, psnr={psnr:.3f} dB")
        detail = f"max {maximum:.3g}; mean {mean:.3g}; PSNR {psnr:.2f} dB"
        if explicitly_tolerated and not exact_enough:
            detail += (
                f"; explicit mean <= {image_mean_tolerance:.3g}, PSNR >= {image_psnr_min:.2f} dB"
            )
        lines.append(f"| `{name}` image | {'pass' if equal else 'fail'} | {detail} |")

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
    snapshot_parser.add_argument(
        "--backend",
        choices=("plotnine", "matplotlib"),
        default="plotnine",
    )
    compare_parser = subparsers.add_parser("compare")
    compare_parser.add_argument("baseline", type=Path)
    compare_parser.add_argument("candidate", type=Path)
    compare_parser.add_argument("--report", type=Path)
    compare_parser.add_argument(
        "--allow-image-difference",
        action="append",
        default=[],
        metavar="NAME",
        help="Allow NAME to use the explicit mean/PSNR tolerances; repeat as needed.",
    )
    compare_parser.add_argument("--image-mean-tolerance", type=float, default=1e-4)
    compare_parser.add_argument("--image-psnr-min", type=float, default=math.inf)
    compare_parser.add_argument("--allow-backend-difference", action="store_true")
    return parser


def main() -> int:
    args = _parser().parse_args()
    if args.command == "snapshot":
        snapshot(args.output, backend=args.backend)
        return 0
    return (
        0
        if compare(
            args.baseline,
            args.candidate,
            args.report,
            allowed_image_differences=set(args.allow_image_difference),
            image_mean_tolerance=args.image_mean_tolerance,
            image_psnr_min=args.image_psnr_min,
            allow_backend_difference=args.allow_backend_difference,
        )
        else 1
    )


if __name__ == "__main__":
    raise SystemExit(main())
