import matplotlib
import pytest

matplotlib.use("Agg")


@pytest.fixture(scope="session")
def adata():
    """Scanpy's bundled PBMC68k subset, which needs no download.

    The example and gallery scripts render their committed figures from
    ``pbmc3k_processed`` instead. This fixture stays on the bundled object so
    the suite runs offline and fast; it gains the same derived ``depth`` column
    those scripts use, so their builders can be smoke-tested here unchanged.
    """
    import numpy as np
    import pandas as pd
    import scanpy as sc

    adata = sc.datasets.pbmc68k_reduced()
    adata.obs["depth"] = pd.Categorical(
        np.where(adata.obs["n_counts"] < adata.obs["n_counts"].median(), "low", "high"),
        categories=["low", "high"],
        ordered=True,
    )
    return adata


@pytest.fixture(scope="session")
def markers(adata):
    candidates = ["CD3D", "NKG7", "CST3", "GNLY", "MS4A1", "FCGR3A", "CD8A"]
    return [g for g in candidates if g in adata.raw.var_names][:4]


@pytest.fixture(scope="session")
def group_key():
    return "bulk_labels"
