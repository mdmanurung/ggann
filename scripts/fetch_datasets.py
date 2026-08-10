"""Download the Scanpy datasets the vignettes and documentation build need.

Run this once from a source checkout while online:

    python scripts/fetch_datasets.py

The files land in ``<repo>/data`` and are reused by every later run. Continuous
integration calls this before the offline vignette and documentation steps, so
those steps never reach the network themselves.
"""

from __future__ import annotations

from pathlib import Path

import scanpy as sc

DATASET_DIR = Path(__file__).resolve().parents[1] / "data"


def main() -> int:
    DATASET_DIR.mkdir(parents=True, exist_ok=True)
    sc.settings.datasetdir = DATASET_DIR
    adata = sc.datasets.pbmc3k_processed()
    print(f"pbmc3k_processed: {adata.n_obs} cells x {adata.n_vars} genes -> {DATASET_DIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
