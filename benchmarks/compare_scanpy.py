#!/usr/bin/env python3
"""Run reproducible, matched ggann-versus-Scanpy plot benchmarks.

The runner keeps correctness and speed claims separate.  Each case validates
the tabular/statistical payload used by both implementations on one immutable
AnnData object before it reports a speed ratio.  Preparation, public plot-object
construction, PNG rendering from a materialized figure, and complete
construct-to-PNG execution are timed independently.  Because ggann and Scanpy
have different deferred-rendering APIs, only preparation and end-to-end ratios
are suitable for release claims; the JSON records that limitation explicitly.
"""

from __future__ import annotations

import argparse
import base64
import gc
import hashlib
import importlib.metadata
import json
import math
import os
import platform
import statistics
import subprocess
import sys
import tempfile
import time
import warnings
from dataclasses import asdict, dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, Callable, Iterable

for _name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_name, "1")
os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "ggann-benchmark-mpl"))
os.environ.setdefault("NUMBA_CACHE_DIR", str(Path(tempfile.gettempdir()) / "ggann-benchmark-numba"))


try:
    from benchmarks.run_benchmarks import (
        CaseSpec,
        _ggann_source_metadata,
        _make_fixture,
        _positive_int_values,
        _rss_bytes,
        _RSSSampler,
        _shape_for,
        _shape_variants,
    )
except ModuleNotFoundError:  # direct ``python benchmarks/compare_scanpy.py``
    from run_benchmarks import (  # type: ignore[no-redef]
        CaseSpec,
        _ggann_source_metadata,
        _make_fixture,
        _positive_int_values,
        _rss_bytes,
        _RSSSampler,
        _shape_for,
        _shape_variants,
    )


_RESULT_PREFIX = "GGANN_SCANPY_RESULT="
_MEMORY_PREFIX = "GGANN_SCANPY_MEMORY="
_FIGSIZE = (6.0, 4.5)
_DPI = 80
_RANK_N_GENES = 2
_SCALING_FIELDS = ("n_obs", "n_vars", "n_genes", "n_groups")

_PRIMARY_WORKLOADS = (
    "embedding_categorical",
    "dotplot",
    "matrixplot",
)
_WORKLOADS = (
    "embedding_categorical",
    "embedding_continuous",
    "embedding_gene",
    "dotplot",
    "matrixplot",
    "violin",
    "stacked_violin",
    "tracksplot",
    "highest_expr_genes",
    "rank_genes_dotplot",
    "rank_genes_matrixplot",
)
_OBS_ONLY_WORKLOADS = {"embedding_categorical", "embedding_continuous"}
_RANK_WORKLOADS = {"rank_genes_dotplot", "rank_genes_matrixplot"}

_STAGE_BOUNDARIES = {
    "preparation": (
        "Library-native extraction and statistical preparation only; no artists. "
        "The payload is validated numerically before ratios are admitted."
    ),
    "construction": (
        "One public plotting call with display disabled. Scanpy creates matplotlib "
        "artists eagerly for embedding, violin, tracksplot, and highest-expression "
        "plots, while ggann returns a deferred plotnine object. Record this stage "
        "for diagnostics, not cross-library speed claims."
    ),
    "render": (
        "PNG save from an already materialized matplotlib Figure at matched size "
        "and DPI. Plot-object construction and artist materialization are excluded."
    ),
    "end_to_end": (
        "Public plotting call, artist materialization, and PNG save at matched size "
        "and DPI. This is the authoritative plotting-speed comparison."
    ),
}


@dataclass(frozen=True)
class ComparisonSpec:
    preset: str
    matrix_format: str
    workload: str
    source: str
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
    ggann_backend: str = "plotnine"

    @property
    def case_id(self) -> str:
        defaults = _shape_for(self.preset, self.matrix_format)
        changed = [
            f"{name}={getattr(self, name)}"
            for name in _SCALING_FIELDS
            if getattr(self, name) != defaults[name]
        ]
        suffix = f"[{','.join(changed)}]" if changed else ""
        backend = (
            f"[ggann_backend={self.ggann_backend}]" if self.ggann_backend != "plotnine" else ""
        )
        return f"{self.preset}{suffix}/{self.matrix_format}/{self.source}/{self.workload}{backend}"

    def fixture_spec(self) -> CaseSpec:
        return CaseSpec(
            preset=self.preset,
            matrix_format=self.matrix_format,
            workload="resolve_x",
            n_obs=self.n_obs,
            n_vars=self.n_vars,
            n_genes=self.n_genes,
            n_groups=self.n_groups,
            density=self.density,
            embedding_dims=self.embedding_dims,
            render_cells=self.render_cells,
            seed=self.seed,
            repeats=self.repeats,
            rss_interval_seconds=self.rss_interval_seconds,
        )


def _selected_workloads(value: str) -> list[str]:
    if value == "primary":
        return list(_PRIMARY_WORKLOADS)
    if value == "all":
        return list(_WORKLOADS)
    selected = [item.strip() for item in value.split(",") if item.strip()]
    if not selected:
        raise ValueError("At least one workload is required.")
    invalid = sorted(set(selected) - set(_WORKLOADS))
    if invalid:
        raise ValueError(f"Unknown workload(s): {', '.join(invalid)}")
    return list(dict.fromkeys(selected))


def _selected_formats(value: str) -> list[str]:
    selected = [item.strip().lower() for item in value.split(",") if item.strip()]
    if selected == ["all"]:
        selected = ["dense", "csr", "csc"]
    if not selected:
        raise ValueError("At least one matrix format is required.")
    invalid = sorted(set(selected) - {"dense", "csr", "csc"})
    if invalid:
        raise ValueError(f"Unknown matrix format(s): {', '.join(invalid)}")
    return list(dict.fromkeys(selected))


def _selected_sources(value: str) -> list[str]:
    selected = [item.strip().lower() for item in value.split(",") if item.strip()]
    if selected == ["all"]:
        selected = ["x", "layer", "raw"]
    if not selected:
        raise ValueError("At least one expression source is required.")
    invalid = sorted(set(selected) - {"x", "layer", "raw"})
    if invalid:
        raise ValueError(f"Unknown source(s): {', '.join(invalid)}")
    return list(dict.fromkeys(selected))


def _source_kwargs(source: str) -> dict[str, Any]:
    if source == "x":
        return {"use_raw": False}
    if source == "layer":
        return {"layer": "counts", "use_raw": False}
    if source == "raw":
        return {"use_raw": True}
    raise ValueError(f"Unknown expression source: {source}")


def _embedding_color(workload: str, genes: list[str]) -> str:
    if workload == "embedding_categorical":
        return "group"
    if workload == "embedding_continuous":
        return "score"
    if workload == "embedding_gene":
        return genes[0]
    raise ValueError(f"Not an embedding workload: {workload}")


def _ordered_unique(values: Iterable[str]) -> list[str]:
    return list(dict.fromkeys(str(value) for value in values))


def _add_fixture_metadata(adata: Any) -> None:
    """Install deterministic colors so neither plotting library mutates ``uns``."""
    from matplotlib import colormaps, colors

    categories = list(adata.obs["group"].cat.categories)
    cmap = colormaps.get_cmap("tab20")
    adata.uns["group_colors"] = [
        colors.to_hex(cmap(index / max(1, len(categories) - 1))) for index in range(len(categories))
    ]


def _ensure_rank_results(adata: Any) -> None:
    import scanpy as sc

    sc.tl.rank_genes_groups(
        adata,
        groupby="group",
        method="t-test",
        use_raw=False,
        n_genes=min(adata.n_vars, max(10, _RANK_N_GENES * adata.obs["group"].nunique())),
        key_added="rank_genes_groups",
    )


def _hash_array(digest: Any, value: Any) -> None:
    import numpy as np
    from scipy import sparse

    if sparse.issparse(value):
        matrix = value
        digest.update(f"sparse:{matrix.format}:{matrix.shape}:{matrix.dtype}".encode())
        for array in (matrix.data, matrix.indices, matrix.indptr):
            contiguous = np.ascontiguousarray(array)
            digest.update(contiguous.view(np.uint8))
        return
    array = np.asarray(value)
    digest.update(f"array:{array.shape}:{array.dtype}".encode())
    if array.dtype.hasobject:
        digest.update(repr(array.tolist()).encode())
    else:
        digest.update(np.ascontiguousarray(array).view(np.uint8))


def _hash_value(digest: Any, value: Any) -> None:
    import numpy as np
    import pandas as pd
    from scipy import sparse

    if isinstance(value, pd.DataFrame):
        digest.update(repr(list(value.columns)).encode())
        digest.update(repr([str(dtype) for dtype in value.dtypes]).encode())
        digest.update(pd.util.hash_pandas_object(value, index=True).to_numpy().tobytes())
    elif isinstance(value, pd.Series):
        digest.update(pd.util.hash_pandas_object(value, index=True).to_numpy().tobytes())
    elif sparse.issparse(value) or isinstance(value, np.ndarray):
        _hash_array(digest, value)
    elif isinstance(value, dict):
        for key in sorted(value, key=str):
            digest.update(f"key:{key!r}".encode())
            _hash_value(digest, value[key])
    elif isinstance(value, (list, tuple)):
        for item in value:
            _hash_value(digest, item)
    elif value is None or isinstance(value, (str, bytes, bool, int, float)):
        digest.update(repr(value).encode())
    else:
        digest.update(repr(value).encode())


def _adata_digest(adata: Any) -> str:
    """Fingerprint all plot-relevant AnnData state without materializing sparsity."""
    digest = hashlib.sha256()
    _hash_array(digest, adata.X)
    _hash_value(digest, adata.obs)
    _hash_value(digest, adata.var)
    # AnnData 0.13 exposes ``None`` in this mapping as an alias for ``.X``.
    # ``.X`` is already hashed above, so only hash named layers here.
    for key in sorted((key for key in adata.layers if key is not None), key=str):
        digest.update(f"layer:{key}".encode())
        _hash_array(digest, adata.layers[key])
    for key in sorted(adata.obsm):
        digest.update(f"obsm:{key}".encode())
        _hash_array(digest, adata.obsm[key])
    if adata.raw is not None:
        digest.update(b"raw")
        _hash_array(digest, adata.raw.X)
        _hash_value(digest, adata.raw.var)
    _hash_value(digest, adata.uns)
    return digest.hexdigest()


def _frame_with_obs_name(frame: Any) -> Any:
    result = frame.copy()
    result.insert(0, "obs_name", result.index.astype(str))
    return result.reset_index(drop=True)


def _ggann_highest_frame(adata: Any, n: int, source: str) -> Any:
    import numpy as np
    import pandas as pd
    from scipy import sparse

    from ggann._expression import resolve_source, source_matrix
    from ggann.qc import _expression_totals, _mean_percentages, _percentage_normalization

    kwargs = _source_kwargs(source)
    kind, layer = resolve_source(adata, kwargs.get("layer"), kwargs.get("use_raw"))
    matrix, var_names = source_matrix(adata, kind, layer)
    totals, elementwise = _expression_totals(matrix)
    denominators, included, _ = _percentage_normalization(totals)
    means = _mean_percentages(matrix, denominators, included, elementwise=elementwise)
    top_indices = np.argsort(-means, kind="stable")[: min(n, len(var_names))]
    selected = matrix[:, top_indices]
    selected = selected.toarray() if sparse.issparse(selected) else np.asarray(selected)
    percentages = np.full(selected.shape, np.nan, dtype=float)
    np.divide(selected, denominators[:, None], out=percentages, where=included[:, None])
    percentages *= 100.0
    return pd.DataFrame(percentages, columns=var_names[top_indices]).melt(
        var_name="gene", value_name="percent"
    )


def _scanpy_highest_frame(adata: Any, n: int, source: str) -> Any:
    import numpy as np
    import pandas as pd
    import scanpy as sc
    from scipy import sparse

    if source == "raw":
        raise ValueError("scanpy.pl.highest_expr_genes has no use_raw parameter")
    layer = "counts" if source == "layer" else None
    normalized = sc.pp.normalize_total(adata, target_sum=100, layer=layer, inplace=False)["X"]
    means = normalized.mean(axis=0).A1 if sparse.issparse(normalized) else normalized.mean(axis=0)
    top_indices = np.argsort(means)[::-1][: min(n, adata.n_vars)]
    selected = normalized[:, top_indices]
    selected = selected.toarray() if sparse.issparse(selected) else np.asarray(selected)
    return pd.DataFrame(selected, columns=adata.var_names[top_indices]).melt(
        var_name="gene", value_name="percent"
    )


def _rank_genes(adata: Any, n_genes: int) -> list[str]:
    import scanpy as sc

    frame = sc.get.rank_genes_groups_df(adata, group=None, key="rank_genes_groups")
    selected = frame.groupby("group", observed=True, sort=False).head(n_genes)
    return _ordered_unique(selected["names"])


def _ggann_preparation(adata: Any, genes: list[str], spec: ComparisonSpec) -> Any:

    from ggann._aggregate import (
        aggregate_expression,
        aggregate_expression_native,
        aggregate_means,
        aggregate_means_native,
        tidy_expression,
    )
    from ggann._resolve import obsm, plain_name, resolve_frame
    from ggann.plots import _native_embedding_frame

    kwargs = _source_kwargs(spec.source)
    workload = spec.workload
    if workload.startswith("embedding_"):
        color = _embedding_color(workload, genes)
        if spec.ggann_backend == "matplotlib":
            frame, x_name, y_name, color_name = _native_embedding_frame(
                adata,
                "X_umap",
                color,
                layer=kwargs.get("layer"),
                use_raw=kwargs.get("use_raw"),
            )
            assert color_name is not None
            coordinate_names = [x_name, y_name]
        else:
            coordinates = [obsm("umap", 0), obsm("umap", 1)]
            frame = resolve_frame(adata, [*coordinates, color], **kwargs)
            color_name = plain_name(adata, color)
            coordinate_names = [plain_name(adata, coordinate) for coordinate in coordinates]
        frame = _frame_with_obs_name(frame)
        return frame.rename(
            columns={
                coordinate_names[0]: "x",
                coordinate_names[1]: "y",
                color_name: "color",
            }
        )[["obs_name", "x", "y", "color"]]
    selected = genes
    if workload in _RANK_WORKLOADS:
        selected = _rank_genes(adata, _RANK_N_GENES)
    if workload in {"dotplot", "rank_genes_dotplot"}:
        aggregate = (
            aggregate_expression_native
            if spec.ggann_backend == "matplotlib" and workload == "dotplot"
            else aggregate_expression
        )
        return aggregate(adata, selected, "group", **kwargs)[
            ["group", "feature", "mean_expression", "fraction"]
        ]
    if workload in {"matrixplot", "rank_genes_matrixplot"}:
        standard_scale = "var" if workload == "rank_genes_matrixplot" else None
        aggregate = (
            aggregate_means_native
            if spec.ggann_backend == "matplotlib" and workload == "matrixplot"
            else aggregate_means
        )
        return aggregate(adata, selected, "group", standard_scale=standard_scale, **kwargs)[
            ["group", "feature", "mean_expression"]
        ]
    if workload in {"violin", "stacked_violin", "tracksplot"}:
        frame = tidy_expression(adata, selected, "group", **kwargs)
        return frame[["obs_name", "group", "feature", "value"]]
    if workload == "highest_expr_genes":
        return _ggann_highest_frame(adata, min(20, adata.n_vars), spec.source)
    raise ValueError(f"Unknown workload: {workload}")


def _scanpy_preparation(adata: Any, genes: list[str], spec: ComparisonSpec) -> Any:
    import pandas as pd
    import scanpy as sc

    kwargs = _source_kwargs(spec.source)
    workload = spec.workload
    if workload.startswith("embedding_"):
        color = _embedding_color(workload, genes)
        coords = adata.obsm["X_umap"][:, :2]
        values = sc.get.obs_df(adata, keys=[color], **kwargs)
        frame = pd.DataFrame(coords, index=adata.obs_names, columns=["x", "y"])
        frame["color"] = values[color]
        return _frame_with_obs_name(frame)[["obs_name", "x", "y", "color"]]

    selected = genes
    if workload in _RANK_WORKLOADS:
        selected = _rank_genes(adata, _RANK_N_GENES)
    if workload == "highest_expr_genes":
        return _scanpy_highest_frame(adata, min(20, adata.n_vars), spec.source)

    wide = sc.get.obs_df(adata, keys=["group", *selected], **kwargs)
    if workload in {"dotplot", "rank_genes_dotplot"}:
        indexed = wide.set_index("group")
        means = indexed.groupby(level=0, observed=True, sort=False).mean()
        fractions = (indexed > 0).groupby(level=0, observed=True, sort=False).mean()
        mean_long = (
            means.rename_axis("group")
            .reset_index()
            .melt(id_vars="group", var_name="feature", value_name="mean_expression")
        )
        fraction_long = (
            fractions.rename_axis("group")
            .reset_index()
            .melt(id_vars="group", var_name="feature", value_name="fraction")
        )
        return mean_long.merge(fraction_long, on=["group", "feature"], validate="one_to_one")
    if workload in {"matrixplot", "rank_genes_matrixplot"}:
        means = wide.set_index("group").groupby(level=0, observed=True, sort=False).mean()
        if workload == "rank_genes_matrixplot":
            means -= means.min(axis=0)
            means = (means / means.max(axis=0)).fillna(0)
        return (
            means.rename_axis("group")
            .reset_index()
            .melt(id_vars="group", var_name="feature", value_name="mean_expression")
        )
    if workload in {"violin", "stacked_violin", "tracksplot"}:
        wide = _frame_with_obs_name(wide)
        return wide.melt(
            id_vars=["obs_name", "group"],
            value_vars=selected,
            var_name="feature",
            value_name="value",
        )
    raise ValueError(f"Unknown workload: {workload}")


def _scanpy_native_preparation(adata: Any, genes: list[str], spec: ComparisonSpec) -> Any:
    """Return Scanpy's native wide/statistical payload without validator reshapes."""
    import numpy as np
    import pandas as pd
    import scanpy as sc
    from scipy import sparse

    kwargs = _source_kwargs(spec.source)
    workload = spec.workload
    if workload.startswith("embedding_"):
        # The embedding path already returns one row per observation, matching
        # Scanpy's scatter input without a canonical-only long reshape.
        return _scanpy_preparation(adata, genes, spec)
    if workload == "highest_expr_genes":
        if spec.source == "raw":
            raise ValueError("scanpy.pl.highest_expr_genes has no use_raw parameter")
        layer = "counts" if spec.source == "layer" else None
        normalized = sc.pp.normalize_total(adata, target_sum=100, layer=layer, inplace=False)["X"]
        means = (
            normalized.mean(axis=0).A1 if sparse.issparse(normalized) else normalized.mean(axis=0)
        )
        top_indices = np.argsort(means)[::-1][: min(20, adata.n_vars)]
        selected = normalized[:, top_indices]
        selected = selected.toarray() if sparse.issparse(selected) else np.asarray(selected)
        return pd.DataFrame(selected, columns=adata.var_names[top_indices])

    selected = genes
    if workload in _RANK_WORKLOADS:
        selected = _rank_genes(adata, _RANK_N_GENES)
    wide = sc.get.obs_df(adata, keys=["group", *selected], **kwargs)
    if workload in {"violin", "stacked_violin", "tracksplot"}:
        return wide
    indexed = wide.set_index("group")
    means = indexed.groupby(level=0, observed=True).mean()
    if workload in {"matrixplot", "rank_genes_matrixplot"}:
        if workload == "rank_genes_matrixplot":
            means -= means.min(axis=0)
            means = (means / means.max(axis=0)).fillna(0)
        return means
    if workload in {"dotplot", "rank_genes_dotplot"}:
        fractions = (indexed > 0).groupby(level=0, observed=True).mean()
        return {"mean_expression": means, "fraction": fractions}
    raise ValueError(f"Unknown workload: {workload}")


def _ggann_native_preparation(adata: Any, genes: list[str], spec: ComparisonSpec) -> Any:
    """Return the compact payload consumed by the explicit Matplotlib backend."""
    if spec.ggann_backend != "matplotlib" or spec.workload not in _PRIMARY_WORKLOADS:
        return _ggann_preparation(adata, genes, spec)

    from ggann._aggregate import grouped_expression_native
    from ggann.plots import _native_embedding_frame

    kwargs = _source_kwargs(spec.source)
    if spec.workload.startswith("embedding_"):
        frame, _, _, _ = _native_embedding_frame(
            adata,
            "X_umap",
            _embedding_color(spec.workload, genes),
            layer=kwargs.get("layer"),
            use_raw=kwargs.get("use_raw"),
        )
        return frame
    means, fractions = grouped_expression_native(
        adata,
        genes,
        "group",
        expression_cutoff=0.0 if spec.workload == "dotplot" else None,
        **kwargs,
    )
    if spec.workload == "dotplot":
        assert fractions is not None
        return {"mean_expression": means, "fraction": fractions}
    return means


def _comparison_columns(workload: str) -> tuple[list[str], list[str], list[str]]:
    if workload.startswith("embedding_"):
        return ["obs_name"], ["color"], ["x", "y"]
    if workload in {"dotplot", "rank_genes_dotplot"}:
        return ["group", "feature"], [], ["mean_expression", "fraction"]
    if workload in {"matrixplot", "rank_genes_matrixplot"}:
        return ["group", "feature"], [], ["mean_expression"]
    if workload in {"violin", "stacked_violin", "tracksplot"}:
        return ["obs_name", "group", "feature"], [], ["value"]
    if workload == "highest_expr_genes":
        return ["gene"], [], ["percent"]
    raise ValueError(f"Unknown workload: {workload}")


def _normalized_strings(series: Any) -> Any:
    return series.astype("string").fillna("<NA>")


def _compare_prepared(left: Any, right: Any, workload: str) -> dict[str, Any]:
    import numpy as np
    from pandas.api.types import is_numeric_dtype

    keys, maybe_numeric, numeric = _comparison_columns(workload)
    columns = [*keys, *maybe_numeric, *numeric]
    missing = {
        "ggann": sorted(set(columns) - set(left.columns)),
        "scanpy": sorted(set(columns) - set(right.columns)),
    }
    issues = []
    if missing["ggann"] or missing["scanpy"]:
        issues.append(f"missing comparison columns: {missing}")
        return {"status": "fail", "issues": issues, "max_abs_diff": None}

    left = left.loc[:, columns].copy()
    right = right.loc[:, columns].copy()
    for column in maybe_numeric:
        left_numeric = is_numeric_dtype(left[column])
        right_numeric = is_numeric_dtype(right[column])
        if left_numeric and right_numeric:
            numeric.append(column)
        else:
            keys.append(column)

    for column in keys:
        left[column] = _normalized_strings(left[column])
        right[column] = _normalized_strings(right[column])

    sort_columns = list(keys)
    if workload == "highest_expr_genes":
        sort_columns.append("percent")
    left = left.sort_values(sort_columns, kind="stable").reset_index(drop=True)
    right = right.sort_values(sort_columns, kind="stable").reset_index(drop=True)
    if left.shape != right.shape:
        issues.append(f"prepared shapes differ: {left.shape} != {right.shape}")
        return {"status": "fail", "issues": issues, "max_abs_diff": None}

    if workload == "highest_expr_genes":
        left_undefined = int(left["percent"].isna().sum())
        right_undefined = int(right["percent"].isna().sum())
        if left_undefined != right_undefined:
            issues.append(
                "undefined percentage rows differ "
                f"(ggann={left_undefined}, scanpy={right_undefined}); "
                "zero-total-cell handling is not semantically matched"
            )

    for column in keys:
        if not _normalized_strings(left[column]).equals(_normalized_strings(right[column])):
            issues.append(f"prepared key column {column!r} differs")

    max_abs_diff = 0.0
    for column in _ordered_unique(numeric):
        left_values = left[column].to_numpy(dtype=float)
        right_values = right[column].to_numpy(dtype=float)
        finite = np.isfinite(left_values) & np.isfinite(right_values)
        if finite.any():
            max_abs_diff = max(
                max_abs_diff,
                float(np.max(np.abs(left_values[finite] - right_values[finite]))),
            )
        if not np.allclose(left_values, right_values, rtol=1e-6, atol=1e-7, equal_nan=True):
            issues.append(f"prepared numeric column {column!r} differs")

    return {
        "status": "pass" if not issues else "fail",
        "issues": issues,
        "ggann_shape": list(left.shape),
        "scanpy_shape": list(right.shape),
        "max_abs_diff": max_abs_diff,
        "rtol": 1e-6,
        "atol": 1e-7,
    }


def _categories(adata: Any) -> list[str]:
    return [str(value) for value in adata.obs["group"].cat.categories]


def _construct(library: str, adata: Any, genes: list[str], spec: ComparisonSpec) -> Any:
    kwargs = _source_kwargs(spec.source)
    workload = spec.workload
    categories = _categories(adata)

    if library == "ggann":
        import ggann as ag

        if workload.startswith("embedding_"):
            return ag.plot_embedding(
                adata,
                "umap",
                color=_embedding_color(workload, genes),
                pointdensity=False,
                downsample=None,
                backend=spec.ggann_backend,
                **kwargs,
            )
        if workload == "dotplot":
            return ag.plot_dotplot(
                adata,
                genes,
                "group",
                categories_order=categories,
                backend=spec.ggann_backend,
                **kwargs,
            )
        if workload == "matrixplot":
            return ag.plot_matrixplot(
                adata,
                genes,
                "group",
                categories_order=categories,
                backend=spec.ggann_backend,
                **kwargs,
            )
        if workload == "violin":
            return ag.plot_violin(
                adata,
                genes,
                "group",
                add_box=True,
                add_points=False,
                downsample=None,
                categories_order=categories,
                **kwargs,
            )
        if workload == "stacked_violin":
            return ag.plot_stacked_violin(
                adata,
                genes,
                "group",
                downsample=None,
                categories_order=categories,
                **kwargs,
            )
        if workload == "tracksplot":
            return ag.plot_tracksplot(adata, genes, "group", categories_order=categories, **kwargs)
        if workload == "highest_expr_genes":
            return ag.plot_highest_expr_genes(adata, n=min(20, adata.n_vars), **kwargs)
        if workload == "rank_genes_dotplot":
            return ag.plot_rank_genes_dotplot(
                adata,
                n_genes=_RANK_N_GENES,
                group_by="group",
                categories_order=categories,
                **kwargs,
            )
        if workload == "rank_genes_matrixplot":
            return ag.plot_rank_genes_matrixplot(
                adata,
                n_genes=_RANK_N_GENES,
                group_by="group",
                standard_scale="var",
                categories_order=categories,
                **kwargs,
            )

    if library == "scanpy":
        import scanpy as sc

        if workload.startswith("embedding_"):
            return sc.pl.embedding(
                adata,
                "umap",
                color=_embedding_color(workload, genes),
                size=4,
                frameon=True,
                show=False,
                return_fig=True,
                **kwargs,
            )
        if workload == "dotplot":
            return sc.pl.dotplot(
                adata,
                genes,
                "group",
                categories_order=categories,
                figsize=_FIGSIZE,
                show=False,
                return_fig=True,
                **kwargs,
            )
        if workload == "matrixplot":
            return sc.pl.matrixplot(
                adata,
                genes,
                "group",
                categories_order=categories,
                figsize=_FIGSIZE,
                show=False,
                return_fig=True,
                **kwargs,
            )
        if workload == "violin":
            return sc.pl.violin(
                adata,
                genes,
                "group",
                order=categories,
                stripplot=False,
                inner="box",
                density_norm="width",
                show=False,
                **kwargs,
            )
        if workload == "stacked_violin":
            return sc.pl.stacked_violin(
                adata,
                genes,
                "group",
                categories_order=categories,
                stripplot=False,
                density_norm="width",
                figsize=_FIGSIZE,
                show=False,
                return_fig=True,
                **kwargs,
            )
        if workload == "tracksplot":
            return sc.pl.tracksplot(
                adata,
                genes,
                "group",
                figsize=_FIGSIZE,
                show=False,
                **kwargs,
            )
        if workload == "highest_expr_genes":
            layer = "counts" if spec.source == "layer" else None
            return sc.pl.highest_expr_genes(
                adata,
                n_top=min(20, adata.n_vars),
                layer=layer,
                show=False,
            )
        rank_kwargs = dict(kwargs)
        if workload == "rank_genes_dotplot":
            return sc.pl.rank_genes_groups_dotplot(
                adata,
                n_genes=_RANK_N_GENES,
                groupby="group",
                key="rank_genes_groups",
                categories_order=categories,
                dendrogram=False,
                figsize=_FIGSIZE,
                show=False,
                return_fig=True,
                **rank_kwargs,
            )
        if workload == "rank_genes_matrixplot":
            return sc.pl.rank_genes_groups_matrixplot(
                adata,
                n_genes=_RANK_N_GENES,
                groupby="group",
                key="rank_genes_groups",
                standard_scale="var",
                categories_order=categories,
                dendrogram=False,
                figsize=_FIGSIZE,
                show=False,
                return_fig=True,
                **rank_kwargs,
            )
    raise ValueError(f"Unknown library/workload pair: {library}/{workload}")


def _figure_from_output(library: str, output: Any) -> Any:
    from matplotlib.figure import Figure

    if library == "ggann":
        return output.draw(show=False)
    if isinstance(output, Figure):
        return output
    if hasattr(output, "make_figure"):
        if getattr(output, "fig", None) is None:
            output.make_figure()
        return output.fig
    figure = getattr(output, "figure", None)
    if figure is not None:
        return figure
    figure = getattr(output, "fig", None)
    if figure is not None:
        return figure
    values = output.values() if isinstance(output, dict) else output
    if isinstance(values, Iterable) and not isinstance(values, (str, bytes)):
        for value in values:
            figure = getattr(value, "figure", None)
            if figure is not None:
                return figure
    raise TypeError(f"Could not recover a matplotlib Figure from {type(output)!r}")


def _close_output(output: Any) -> None:
    from matplotlib import pyplot as plt

    figures = []
    for attribute in ("figure", "fig"):
        figure = getattr(output, attribute, None)
        if figure is not None:
            figures.append(figure)
    if isinstance(output, dict):
        figures.extend(value.figure for value in output.values() if hasattr(value, "figure"))
    elif isinstance(output, (list, tuple)):
        figures.extend(value.figure for value in output if hasattr(value, "figure"))
    for figure in dict.fromkeys(figures):
        plt.close(figure)


def _output_description(output: Any) -> dict[str, Any]:
    import pandas as pd

    description: dict[str, Any] = {"kind": f"{type(output).__module__}.{type(output).__name__}"}
    if isinstance(output, pd.DataFrame):
        description["shape"] = list(output.shape)
        description["bytes"] = int(output.memory_usage(index=True, deep=True).sum())
    elif isinstance(output, dict):
        description["keys"] = sorted(str(key) for key in output)
    elif hasattr(output, "data") and isinstance(output.data, pd.DataFrame):
        description["plot_data_shape"] = list(output.data.shape)
    return description


def _measure_call(
    function: Callable[[], Any],
    rss_interval_seconds: float,
    *,
    cleanup: Callable[[Any], None] | None = None,
) -> dict[str, Any]:
    gc.collect()
    baseline_rss = _rss_bytes()
    with _RSSSampler(rss_interval_seconds) as sampler:
        start = time.perf_counter_ns()
        output = function()
        duration_ns = time.perf_counter_ns() - start
    output_rss = _rss_bytes()
    peak_rss = max(sampler.peak_bytes, output_rss)
    description = _output_description(output)
    if cleanup is not None:
        cleanup(output)
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


def _summarize_samples(samples: list[dict[str, Any]]) -> dict[str, Any]:
    durations = [sample["duration_seconds"] for sample in samples]
    return {
        "samples": samples,
        "median_duration_seconds": statistics.median(durations),
        "min_duration_seconds": min(durations),
        "max_duration_seconds": max(durations),
        "max_peak_rss_delta_bytes": max(sample["peak_rss_delta_bytes"] for sample in samples),
        "max_retained_after_gc_bytes": max(sample["retained_after_gc_bytes"] for sample in samples),
    }


def _measure_pair(
    functions: dict[str, Callable[[], Any]],
    spec: ComparisonSpec,
    *,
    cleanup: Callable[[Any], None] | None = None,
) -> dict[str, Any]:
    libraries = ("ggann", "scanpy")
    cold = {
        library: _measure_call(functions[library], spec.rss_interval_seconds, cleanup=cleanup)
        for library in libraries
    }
    samples = {library: [] for library in libraries}
    execution_order = []
    for repeat in range(spec.repeats):
        order = libraries if repeat % 2 == 0 else tuple(reversed(libraries))
        execution_order.append(list(order))
        for library in order:
            samples[library].append(
                _measure_call(functions[library], spec.rss_interval_seconds, cleanup=cleanup)
            )
    summaries = {
        library: {
            "cold": cold[library],
            "repeated": _summarize_samples(samples[library]),
        }
        for library in libraries
    }
    ggann_median = summaries["ggann"]["repeated"]["median_duration_seconds"]
    scanpy_median = summaries["scanpy"]["repeated"]["median_duration_seconds"]
    return {
        "libraries": summaries,
        "execution_order": execution_order,
        "speedup_scanpy_over_ggann": (scanpy_median / ggann_median if ggann_median else math.inf),
    }


def _prepare_render_function(
    library: str,
    adata: Any,
    genes: list[str],
    spec: ComparisonSpec,
    output_path: Path,
) -> tuple[Callable[[], Any], Callable[[], None]]:
    output = _construct(library, adata, genes, spec)
    figure = _figure_from_output(library, output)
    figure.set_size_inches(*_FIGSIZE, forward=True)

    def render() -> dict[str, int]:
        figure.savefig(output_path, format="png", dpi=_DPI)
        return {"png_bytes": output_path.stat().st_size, "axes": len(figure.axes)}

    def cleanup() -> None:
        from matplotlib import pyplot as plt

        plt.close(figure)
        _close_output(output)

    return render, cleanup


def _end_to_end_function(
    library: str,
    adata: Any,
    genes: list[str],
    spec: ComparisonSpec,
    output_path: Path,
) -> Callable[[], Any]:
    def execute() -> dict[str, int]:
        from matplotlib import pyplot as plt

        output = _construct(library, adata, genes, spec)
        figure = _figure_from_output(library, output)
        figure.set_size_inches(*_FIGSIZE, forward=True)
        figure.savefig(output_path, format="png", dpi=_DPI)
        result = {"png_bytes": output_path.stat().st_size, "axes": len(figure.axes)}
        plt.close(figure)
        _close_output(output)
        return result

    return execute


def _stage_claim_status(stage: str, validation: dict[str, Any]) -> str:
    if validation["status"] != "pass":
        return "not_comparable"
    if stage in {"preparation", "end_to_end"}:
        return "comparable"
    return "diagnostic_only"


def _execute_case(spec: ComparisonSpec) -> dict[str, Any]:
    import matplotlib

    matplotlib.use("Agg", force=True)
    adata, genes, input_bytes = _make_fixture(spec.fixture_spec())
    _add_fixture_metadata(adata)
    if spec.workload in _RANK_WORKLOADS:
        _ensure_rank_results(adata)
    initial_digest = _adata_digest(adata)

    prepared = {
        "ggann": _ggann_preparation(adata, genes, spec),
        "scanpy": _scanpy_preparation(adata, genes, spec),
    }
    validation = _compare_prepared(prepared["ggann"], prepared["scanpy"], spec.workload)

    with tempfile.TemporaryDirectory(prefix="ggann-scanpy-") as temporary:
        temporary_path = Path(temporary)
        stages: dict[str, Any] = {}

        stages["preparation"] = _measure_pair(
            {
                "ggann": lambda: _ggann_native_preparation(adata, genes, spec),
                "scanpy": lambda: _scanpy_native_preparation(adata, genes, spec),
            },
            spec,
        )
        stages["construction"] = _measure_pair(
            {
                library: (lambda library=library: _construct(library, adata, genes, spec))
                for library in ("ggann", "scanpy")
            },
            spec,
            cleanup=_close_output,
        )

        render_functions = {}
        render_cleanups = []
        for library in ("ggann", "scanpy"):
            function, cleanup = _prepare_render_function(
                library,
                adata,
                genes,
                spec,
                temporary_path / f"{library}-render.png",
            )
            render_functions[library] = function
            render_cleanups.append(cleanup)
        try:
            stages["render"] = _measure_pair(render_functions, spec)
        finally:
            for cleanup in render_cleanups:
                cleanup()

        stages["end_to_end"] = _measure_pair(
            {
                library: _end_to_end_function(
                    library,
                    adata,
                    genes,
                    spec,
                    temporary_path / f"{library}-end-to-end.png",
                )
                for library in ("ggann", "scanpy")
            },
            spec,
        )

    final_digest = _adata_digest(adata)
    immutable = initial_digest == final_digest
    if not immutable:
        validation["status"] = "fail"
        validation["issues"].append("AnnData fingerprint changed during benchmark")
    for stage, result in stages.items():
        result["claim_status"] = _stage_claim_status(stage, validation)

    return {
        "case_id": spec.case_id,
        "parameters": asdict(spec),
        "selected_genes": genes,
        "input_bytes": input_bytes,
        "input_fingerprint_before": initial_digest,
        "input_fingerprint_after": final_digest,
        "input_immutable": immutable,
        "comparability": validation,
        "stage_boundaries": _STAGE_BOUNDARIES,
        "stages": stages,
    }


def _execute_isolated_memory_probe(
    spec: ComparisonSpec, library: str, stage: str
) -> dict[str, Any]:
    """Measure one first call in a fresh process for interpretable RSS deltas."""
    import matplotlib

    matplotlib.use("Agg", force=True)
    adata, genes, input_bytes = _make_fixture(spec.fixture_spec())
    _add_fixture_metadata(adata)
    if spec.workload in _RANK_WORKLOADS:
        _ensure_rank_results(adata)
    # Package-import cold starts are measured by ``_cold_import_probe``. Keep
    # imports outside this RSS baseline so the memory ratio describes the plot
    # workload rather than the very different dependency trees.
    if library == "ggann":
        import ggann  # noqa: F401
    elif library == "scanpy":
        import scanpy  # noqa: F401
    else:
        raise ValueError(f"Unknown library: {library}")
    gc.collect()
    before = _adata_digest(adata)
    if stage == "preparation":
        function = (
            (lambda: _ggann_native_preparation(adata, genes, spec))
            if library == "ggann"
            else (lambda: _scanpy_native_preparation(adata, genes, spec))
        )
    elif stage == "end_to_end":
        temporary = tempfile.TemporaryDirectory(prefix="ggann-scanpy-memory-")
        function = _end_to_end_function(
            library,
            adata,
            genes,
            spec,
            Path(temporary.name) / f"{library}.png",
        )
    else:
        raise ValueError(f"Unsupported isolated-memory stage: {stage}")
    try:
        sample = _measure_call(function, spec.rss_interval_seconds)
    finally:
        if stage == "end_to_end":
            temporary.cleanup()
    after = _adata_digest(adata)
    return {
        "library": library,
        "stage": stage,
        "sample": sample,
        "input_bytes": input_bytes,
        "input_fingerprint_before": before,
        "input_fingerprint_after": after,
        "input_immutable": before == after,
        "process_scope": (
            "fresh child; imports and fixture creation precede the RSS baseline; "
            "the measured call is the first plotting/preparation call"
        ),
    }


def _encoded_case(spec: ComparisonSpec) -> str:
    return base64.urlsafe_b64encode(
        json.dumps(asdict(spec), separators=(",", ":")).encode()
    ).decode()


def _decoded_case(payload: str) -> ComparisonSpec:
    return ComparisonSpec(**json.loads(base64.urlsafe_b64decode(payload.encode())))


def _encoded_memory_probe(spec: ComparisonSpec, library: str, stage: str) -> str:
    payload = {"spec": asdict(spec), "library": library, "stage": stage}
    return base64.urlsafe_b64encode(json.dumps(payload, separators=(",", ":")).encode()).decode()


def _run_case_process(spec: ComparisonSpec, timeout_seconds: float) -> dict[str, Any]:
    environment = os.environ.copy()
    environment.update(
        {
            "PYTHONHASHSEED": "0",
            "MPLBACKEND": "Agg",
            "MPLCONFIGDIR": os.environ["MPLCONFIGDIR"],
            "NUMBA_CACHE_DIR": os.environ["NUMBA_CACHE_DIR"],
        }
    )
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
            f"Scanpy comparison case {spec.case_id} failed with exit code "
            f"{completed.returncode}:\n{completed.stderr.strip()}"
        )
    lines = [
        line.removeprefix(_RESULT_PREFIX)
        for line in completed.stdout.splitlines()
        if line.startswith(_RESULT_PREFIX)
    ]
    if len(lines) != 1:
        raise RuntimeError(
            f"Case {spec.case_id} returned no parseable result.\n"
            f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
        )
    result = json.loads(lines[0])
    diagnostics = [line for line in completed.stderr.splitlines() if line.strip()]
    if diagnostics:
        result["diagnostics"] = diagnostics
    return result


def _run_isolated_memory_process(
    spec: ComparisonSpec,
    library: str,
    stage: str,
    timeout_seconds: float,
) -> dict[str, Any]:
    environment = os.environ.copy()
    environment.update(
        {
            "PYTHONHASHSEED": "0",
            "MPLBACKEND": "Agg",
            "MPLCONFIGDIR": os.environ["MPLCONFIGDIR"],
            "NUMBA_CACHE_DIR": os.environ["NUMBA_CACHE_DIR"],
        }
    )
    completed = subprocess.run(
        [
            sys.executable,
            str(Path(__file__).resolve()),
            "--internal-memory",
            _encoded_memory_probe(spec, library, stage),
        ],
        check=False,
        capture_output=True,
        text=True,
        env=environment,
        timeout=timeout_seconds,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            f"Isolated memory probe {spec.case_id}/{library}/{stage} failed "
            f"with exit code {completed.returncode}:\n{completed.stderr.strip()}"
        )
    lines = [
        line.removeprefix(_MEMORY_PREFIX)
        for line in completed.stdout.splitlines()
        if line.startswith(_MEMORY_PREFIX)
    ]
    if len(lines) != 1:
        raise RuntimeError(
            f"Isolated memory probe {spec.case_id}/{library}/{stage} returned "
            f"no parseable result.\nstdout:\n{completed.stdout}\n"
            f"stderr:\n{completed.stderr}"
        )
    result = json.loads(lines[0])
    diagnostics = [line for line in completed.stderr.splitlines() if line.strip()]
    if diagnostics:
        result["diagnostics"] = diagnostics
    return result


def _summarize_isolated_memory_probes(probes: list[dict[str, Any]]) -> dict[str, Any]:
    """Retain raw child results and summarize their independently measured deltas."""
    if not probes:
        raise ValueError("At least one isolated-memory probe is required.")

    def metric(name: str) -> dict[str, Any]:
        values = [probe["sample"][name] for probe in probes]
        return {
            "samples": values,
            "median": statistics.median(values),
            "min": min(values),
            "max": max(values),
        }

    first = probes[0]
    return {
        "library": first["library"],
        "stage": first["stage"],
        "repeat_count": len(probes),
        "probes": probes,
        "summary": {
            name: metric(name)
            for name in (
                "peak_rss_delta_bytes",
                "retained_with_output_bytes",
                "retained_after_gc_bytes",
            )
        },
        "input_immutable": all(probe["input_immutable"] for probe in probes),
        "process_scope": (
            "independent fresh children; imports and fixture creation precede each RSS "
            "baseline; summaries use the median and retain every raw range/sample"
        ),
    }


def _isolated_memory_metric(record: dict[str, Any], name: str) -> float:
    """Read the repeated-child median or a legacy single-child sample."""
    if "summary" in record:
        return float(record["summary"][name]["median"])
    return float(record["sample"][name])


def _versions() -> dict[str, str | None]:
    packages = (
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
        "seaborn",
        "pillow",
    )
    versions = {}
    for package in packages:
        try:
            versions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            versions[package] = None
    return versions


def _cpu_model() -> str | None:
    cpuinfo = Path("/proc/cpuinfo")
    if cpuinfo.exists():
        for line in cpuinfo.read_text().splitlines():
            if line.lower().startswith("model name"):
                return line.split(":", 1)[1].strip()
    return platform.processor() or None


def _git_metadata() -> dict[str, Any]:
    def git(*arguments: str) -> str | None:
        completed = subprocess.run(["git", *arguments], check=False, capture_output=True, text=True)
        return completed.stdout.strip() if completed.returncode == 0 else None

    status = git("status", "--porcelain")
    return {
        "revision": git("rev-parse", "HEAD"),
        "dirty": bool(status) if status is not None else None,
    }


def _cold_import_probe(package: str, timeout_seconds: float) -> dict[str, Any]:
    environment = os.environ.copy()
    environment.update(
        {
            "PYTHONHASHSEED": "0",
            "MPLBACKEND": "Agg",
            "MPLCONFIGDIR": os.environ["MPLCONFIGDIR"],
            "NUMBA_CACHE_DIR": os.environ["NUMBA_CACHE_DIR"],
        }
    )
    start = time.perf_counter()
    completed = subprocess.run(
        [sys.executable, "-c", f"import {package}"],
        check=False,
        capture_output=True,
        text=True,
        env=environment,
        timeout=timeout_seconds,
    )
    return {
        "wall_seconds": time.perf_counter() - start,
        "exit_code": completed.returncode,
        "stderr": completed.stderr.strip() or None,
        "scope": "fresh Python process plus package import; fixture and plotting excluded",
    }


def _geometric_mean(values: list[float]) -> float | None:
    positive = [value for value in values if value > 0 and math.isfinite(value)]
    if len(positive) != len(values) or not positive:
        return None
    return math.exp(sum(math.log(value) for value in positive) / len(positive))


def _evaluate_release_gates(results: list[dict[str, Any]]) -> dict[str, Any]:
    required = set(_PRIMARY_WORKLOADS)
    large_sparse = [
        result
        for result in results
        if result["parameters"]["preset"] == "extended"
        and result["parameters"]["matrix_format"] in {"csr", "csc"}
        and result["parameters"]["source"] == "x"
        and result["parameters"]["workload"] in required
    ]
    by_workload = {result["parameters"]["workload"]: result for result in large_sparse}
    missing = sorted(required - set(by_workload))
    if missing:
        return {
            "status": "not_evaluated",
            "reason": "missing extended sparse primary cases",
            "missing_workloads": missing,
            "baseline_regression_gate": "not_evaluated_by_this_cross-library_runner",
        }

    checks = []
    end_to_end_speedups = []
    for workload in _PRIMARY_WORKLOADS:
        result = by_workload[workload]
        comparable = result["comparability"]["status"] == "pass"
        preparation = result["stages"]["preparation"]
        end_to_end = result["stages"]["end_to_end"]
        prep_speedup = preparation["speedup_scanpy_over_ggann"]
        e2e_speedup = end_to_end["speedup_scanpy_over_ggann"]
        isolated = result.get("isolated_memory", {}).get("end_to_end", {})
        if set(isolated) == {"ggann", "scanpy"}:
            ggann_peak = _isolated_memory_metric(isolated["ggann"], "peak_rss_delta_bytes")
            scanpy_peak = _isolated_memory_metric(isolated["scanpy"], "peak_rss_delta_bytes")
            repeats = isolated["ggann"].get("repeat_count", 1)
            memory_source = f"fresh_child_end_to_end_median_of_{repeats}"
        else:
            ggann_peak = end_to_end["libraries"]["ggann"]["repeated"]["max_peak_rss_delta_bytes"]
            scanpy_peak = end_to_end["libraries"]["scanpy"]["repeated"]["max_peak_rss_delta_bytes"]
            memory_source = "in_case_sampled_rss"
        memory_ratio = ggann_peak / scanpy_peak if scanpy_peak else None
        checks.append(
            {
                "workload": workload,
                "comparable": comparable,
                "preparation_speedup": prep_speedup,
                "preparation_at_least_2x": comparable and prep_speedup >= 2.0,
                "end_to_end_speedup": e2e_speedup,
                "end_to_end_no_more_than_10pct_slower": comparable and e2e_speedup >= (1 / 1.10),
                "peak_memory_ggann_over_scanpy": memory_ratio,
                "peak_memory_source": memory_source,
                "peak_memory_no_more_than_10pct_worse": comparable
                and memory_ratio is not None
                and memory_ratio <= 1.10,
            }
        )
        if comparable:
            end_to_end_speedups.append(e2e_speedup)
    geometric_speedup = _geometric_mean(end_to_end_speedups)
    cross_library_pass = (
        all(
            check["preparation_at_least_2x"]
            and check["end_to_end_no_more_than_10pct_slower"]
            and check["peak_memory_no_more_than_10pct_worse"]
            for check in checks
        )
        and geometric_speedup is not None
        and geometric_speedup > 1.0
    )
    return {
        "status": "pass" if cross_library_pass else "fail",
        "checks": checks,
        "end_to_end_geometric_speedup": geometric_speedup,
        "end_to_end_geometric_mean_faster": geometric_speedup is not None
        and geometric_speedup > 1.0,
        "baseline_regression_gate": "not_evaluated_by_this_cross-library_runner",
        "note": (
            "A cross-library pass is not a complete release-performance pass until "
            "the separate ggann baseline regression gate is also evaluated."
        ),
    }


def _seconds(value: float) -> str:
    if value < 0.001:
        return f"{value * 1_000_000:.0f} us"
    if value < 1:
        return f"{value * 1_000:.1f} ms"
    return f"{value:.3f} s"


def _markdown_report(document: dict[str, Any]) -> str:
    lines = [
        "# ggann versus Scanpy benchmark",
        "",
        "Only rows with `pass` comparability support speed claims. Construction is ",
        "diagnostic because the two public APIs defer different amounts of rendering.",
        "",
        "| Case | Comparable | ggann prep | Scanpy prep | Prep speedup | ggann end-to-end | Scanpy end-to-end | End-to-end speedup |",
        "|---|:---:|---:|---:|---:|---:|---:|---:|",
    ]
    for result in document["results"]:
        prep = result["stages"]["preparation"]
        end = result["stages"]["end_to_end"]
        gp = prep["libraries"]["ggann"]["repeated"]["median_duration_seconds"]
        sp = prep["libraries"]["scanpy"]["repeated"]["median_duration_seconds"]
        ge = end["libraries"]["ggann"]["repeated"]["median_duration_seconds"]
        se = end["libraries"]["scanpy"]["repeated"]["median_duration_seconds"]
        lines.append(
            f"| `{result['case_id']}` | {result['comparability']['status']} | "
            f"{_seconds(gp)} | {_seconds(sp)} | {prep['speedup_scanpy_over_ggann']:.2f}x | "
            f"{_seconds(ge)} | {_seconds(se)} | {end['speedup_scanpy_over_ggann']:.2f}x |"
        )
    lines.extend(
        [
            "",
            "Speedup is Scanpy time divided by ggann time; values above 1 favor ggann.",
            "Raw timing and RSS samples, provenance, input fingerprints, and gate details ",
            "are retained in the JSON document.",
            "",
        ]
    )
    return "\n".join(lines)


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preset", choices=("smoke", "standard", "extended"), default="smoke")
    parser.add_argument("--formats", default="csr", help="dense, csr, csc, or all")
    parser.add_argument(
        "--workloads", default="primary", help="primary, all, or a comma-separated list"
    )
    parser.add_argument("--sources", default="x", help="x, layer, raw, or all")
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument(
        "--ggann-backend",
        choices=("plotnine", "matplotlib"),
        default="plotnine",
        help="ggann rendering/preparation path; Scanpy is unchanged.",
    )
    parser.add_argument("--seed", type=int, default=20_260_809)
    parser.add_argument("--rss-interval-ms", type=float, default=1.0)
    parser.add_argument("--timeout-seconds", type=float, default=1_800.0)
    parser.add_argument("--label", default="ggann-vs-scanpy")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--report", type=Path)
    parser.add_argument("--include-cold-start", action="store_true")
    parser.add_argument(
        "--isolated-memory-stages",
        default="",
        help=(
            "Comma-separated preparation,end_to_end. Each library/stage runs once "
            "in a fresh child; use this for release RSS evidence."
        ),
    )
    parser.add_argument(
        "--isolated-memory-repeats",
        type=int,
        default=1,
        help="Independent fresh children per library and selected memory stage.",
    )
    for field in _SCALING_FIELDS:
        parser.add_argument(
            "--" + field.replace("_", "-"),
            type=_positive_int_values,
            help="Override one preset dimension; one dimension may be a comma-separated sweep.",
        )
    parser.add_argument("--list-workloads", action="store_true")
    parser.add_argument("--internal-case", help=argparse.SUPPRESS)
    parser.add_argument("--internal-memory", help=argparse.SUPPRESS)
    return parser


def _internal_main(payload: str) -> int:
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", FutureWarning)
            result = _execute_case(_decoded_case(payload))
    except Exception:
        import traceback

        traceback.print_exc()
        return 1
    print(_RESULT_PREFIX + json.dumps(result, separators=(",", ":")))
    return 0


def _internal_memory_main(payload: str) -> int:
    try:
        decoded = json.loads(base64.urlsafe_b64decode(payload.encode()))
        spec = ComparisonSpec(**decoded["spec"])
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", FutureWarning)
            result = _execute_isolated_memory_probe(spec, decoded["library"], decoded["stage"])
    except Exception:
        import traceback

        traceback.print_exc()
        return 1
    print(_MEMORY_PREFIX + json.dumps(result, separators=(",", ":")))
    return 0


def _case_supported(workload: str, source: str) -> tuple[bool, str | None]:
    if workload in _OBS_ONLY_WORKLOADS and source != "x":
        return False, "observation-colored embedding is expression-source independent"
    if workload == "highest_expr_genes" and source == "raw":
        return False, "scanpy.pl.highest_expr_genes has no use_raw parameter"
    return True, None


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    if args.internal_case:
        return _internal_main(args.internal_case)
    if args.internal_memory:
        return _internal_memory_main(args.internal_memory)
    if args.list_workloads:
        for workload in _WORKLOADS:
            print(workload)
        return 0
    if args.repeats < 1:
        parser.error("--repeats must be at least 1")
    if args.rss_interval_ms <= 0:
        parser.error("--rss-interval-ms must be positive")
    if args.isolated_memory_repeats < 1:
        parser.error("--isolated-memory-repeats must be at least 1")
    if args.output is None:
        parser.error("--output is required")

    try:
        formats = _selected_formats(args.formats)
        workloads = _selected_workloads(args.workloads)
        sources = _selected_sources(args.sources)
        isolated_memory_stages = [
            item.strip() for item in args.isolated_memory_stages.split(",") if item.strip()
        ]
        invalid_memory_stages = sorted(set(isolated_memory_stages) - {"preparation", "end_to_end"})
        if invalid_memory_stages:
            raise ValueError(
                "Unknown isolated-memory stage(s): " + ", ".join(invalid_memory_stages)
            )
        isolated_memory_stages = list(dict.fromkeys(isolated_memory_stages))
        overrides = {field: getattr(args, field) for field in _SCALING_FIELDS}
        shapes = {
            matrix_format: _shape_variants(args.preset, matrix_format, overrides)
            for matrix_format in formats
        }
    except ValueError as error:
        parser.error(str(error))

    cases = []
    skipped = []
    for matrix_format in formats:
        for shape in shapes[matrix_format]:
            for source in sources:
                for workload in workloads:
                    supported, reason = _case_supported(workload, source)
                    if not supported:
                        skipped.append(
                            {
                                "matrix_format": matrix_format,
                                "source": source,
                                "workload": workload,
                                "reason": reason,
                            }
                        )
                        continue
                    cases.append(
                        ComparisonSpec(
                            preset=args.preset,
                            matrix_format=matrix_format,
                            workload=workload,
                            source=source,
                            seed=args.seed,
                            repeats=args.repeats,
                            rss_interval_seconds=args.rss_interval_ms / 1_000,
                            ggann_backend=args.ggann_backend,
                            **shape,
                        )
                    )

    package_versions = _versions()
    ggann_source = _ggann_source_metadata()
    cold_imports = None
    if args.include_cold_start:
        cold_imports = {
            package: _cold_import_probe(package, args.timeout_seconds)
            for package in ("ggann", "scanpy")
        }

    results = []
    for index, spec in enumerate(cases, start=1):
        print(f"[{index}/{len(cases)}] {spec.case_id}", flush=True)
        result = _run_case_process(spec, args.timeout_seconds)
        if isolated_memory_stages:
            result["isolated_memory"] = {}
            for stage in isolated_memory_stages:
                probes = {library: [] for library in ("ggann", "scanpy")}
                for repeat in range(args.isolated_memory_repeats):
                    order = ("ggann", "scanpy") if repeat % 2 == 0 else ("scanpy", "ggann")
                    for library in order:
                        probes[library].append(
                            _run_isolated_memory_process(
                                spec,
                                library,
                                stage,
                                args.timeout_seconds,
                            )
                        )
                result["isolated_memory"][stage] = {
                    library: _summarize_isolated_memory_probes(probes[library])
                    for library in ("ggann", "scanpy")
                }
        results.append(result)

    document = {
        "schema_version": 1,
        "benchmark_kind": "ggann_scanpy_matched",
        "metadata": {
            "label": args.label,
            "created_at": datetime.now(UTC).isoformat(),
            "python": sys.version,
            "platform": platform.platform(),
            "cpu_model": _cpu_model(),
            "logical_cpus": os.cpu_count(),
            "packages": package_versions,
            "ggann_source": ggann_source,
            "git": _git_metadata(),
            "thread_settings": {
                name: os.environ.get(name)
                for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS")
            },
            "figure": {
                "width_inches": _FIGSIZE[0],
                "height_inches": _FIGSIZE[1],
                "dpi": _DPI,
                "format": "png",
                "backend": "Agg",
            },
            "preset": args.preset,
            "formats": formats,
            "workloads": workloads,
            "sources": sources,
            "shape_overrides": {
                name: value for name, value in overrides.items() if value is not None
            },
            "repeats": args.repeats,
            "ggann_backend": args.ggann_backend,
            "seed": args.seed,
            "rss_interval_ms": args.rss_interval_ms,
            "cold_imports": cold_imports,
            "isolated_memory_stages": isolated_memory_stages,
            "isolated_memory_repeats": args.isolated_memory_repeats,
            "skipped_cases": skipped,
            "claim_policy": (
                "Only preparation and end_to_end stages with comparability=pass may support speed claims. "
                "Construction and render-only stages are diagnostic."
            ),
        },
        "results": results,
        "release_gates": _evaluate_release_gates(results),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(document, indent=2) + "\n")
    if args.report is not None:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(_markdown_report(document))
    print(f"Wrote {len(results)} matched cases to {args.output}")
    return 0 if all(result["comparability"]["status"] == "pass" for result in results) else 2


if __name__ == "__main__":
    raise SystemExit(main())
