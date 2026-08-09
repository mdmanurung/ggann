from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

from benchmarks.render_scanpy_vignette import render_document

ROOT = Path(__file__).resolve().parents[1]
VIGNETTES = sorted((ROOT / "examples" / "vignettes").glob("[0-9]*.py"))


@pytest.mark.parametrize("script", VIGNETTES, ids=lambda path: path.stem)
def test_offline_vignette_executes(script: Path) -> None:
    env = os.environ.copy()
    env.update(
        {
            "GGANN_DOCS_OFFLINE": "1",
            "MPLBACKEND": "Agg",
            "PYTHONHASHSEED": "0",
            "PYTHONPATH": os.pathsep.join([str(ROOT / "src"), env.get("PYTHONPATH", "")]).rstrip(
                os.pathsep
            ),
        }
    )
    result = subprocess.run(
        [sys.executable, str(script)],
        cwd=ROOT,
        env=env,
        text=True,
        capture_output=True,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr


def test_scanpy_vignette_table_matches_recorded_json() -> None:
    document = json.loads((ROOT / "benchmarks/results/scanpy-extended-csr.json").read_text())
    expected = render_document(document)
    actual = (ROOT / "docs/_includes/scanpy-extended-csr.md").read_text()
    assert actual == expected
