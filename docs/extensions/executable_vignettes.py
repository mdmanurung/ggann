"""Execute the deterministic vignette scripts as part of every docs build."""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

from sphinx.errors import ExtensionError
from sphinx.util import logging

LOGGER = logging.getLogger(__name__)
ROOT = Path(__file__).resolve().parents[2]
VIGNETTES = tuple(sorted((ROOT / "examples" / "vignettes").glob("[0-9]*.py")))


def _execute_vignettes(app) -> None:
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
    for path in VIGNETTES:
        LOGGER.info("executing offline vignette: %s", path.name)
        result = subprocess.run(
            [sys.executable, str(path)],
            cwd=ROOT,
            env=env,
            text=True,
            capture_output=True,
            check=False,
        )
        if result.returncode:
            details = "\n".join(part for part in (result.stdout, result.stderr) if part)
            raise ExtensionError(f"Vignette {path.name} failed:\n{details}")


def setup(app):
    app.connect("builder-inited", _execute_vignettes)
    return {
        "version": "1.0",
        "parallel_read_safe": False,
        "parallel_write_safe": True,
    }
