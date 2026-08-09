"""Run Pyright against ggann with the active Python environment on its import path."""

from __future__ import annotations

import json
import site
import subprocess
import sys
import tempfile
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main() -> int:
    """Create a portable transient config and return Pyright's exit status."""
    project = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    config = dict(project["tool"]["pyright"])
    config["include"] = ["../src/ggann"]
    config["extraPaths"] = sorted({str(Path(path).resolve()) for path in site.getsitepackages()})

    with tempfile.TemporaryDirectory(prefix=".pyright-", dir=ROOT) as directory:
        config_path = Path(directory) / "pyrightconfig.json"
        config_path.write_text(json.dumps(config, indent=2), encoding="utf-8")
        return subprocess.call(
            [sys.executable, "-m", "pyright", "--project", str(config_path)],
            cwd=ROOT,
        )


if __name__ == "__main__":
    raise SystemExit(main())
