from __future__ import annotations

import re
import tomllib
from pathlib import Path

from packaging.requirements import Requirement

ROOT = Path(__file__).parents[1]


def _project_metadata() -> dict:
    with (ROOT / "pyproject.toml").open("rb") as handle:
        return tomllib.load(handle)["project"]


def _requirement(name: str) -> Requirement:
    requirements = {
        requirement.name: requirement
        for requirement in map(Requirement, _project_metadata()["dependencies"])
    }
    return requirements[name]


def test_runtime_dependencies_are_index_installable_version_ranges():
    for raw in _project_metadata()["dependencies"]:
        requirement = Requirement(raw)
        assert requirement.url is None, f"direct URL dependency is not releasable: {raw}"
        assert requirement.specifier, f"runtime dependency is unbounded: {raw}"


def test_annplyr_and_anndata_ranges_match_the_v03_contract():
    annplyr = _requirement("annplyr").specifier
    assert annplyr.contains("0.3.0")
    assert not annplyr.contains("0.2.0")
    assert not annplyr.contains("0.4.0")

    anndata = _requirement("anndata").specifier
    assert anndata.contains("0.12.0")
    assert not anndata.contains("0.11.0")


def test_plotnine_extra_range_covers_the_required_public_api():
    versions = _requirement("plotnine-extra").specifier
    assert versions.contains("0.3.1")
    assert not versions.contains("0.3.0")
    assert not versions.contains("0.4.0")


def test_release_support_files_are_present():
    for filename in (
        "CHANGELOG.md",
        "CITATION.cff",
        "CONTRIBUTING.md",
        "LICENSE",
        "README.md",
        "SECURITY.md",
    ):
        path = ROOT / filename
        assert path.is_file() and path.stat().st_size > 0, filename


def test_citation_version_matches_package_version():
    version = _project_metadata()["version"]
    citation = (ROOT / "CITATION.cff").read_text(encoding="utf-8")
    match = re.search(r"^version:\s*([^\s]+)\s*$", citation, flags=re.MULTILINE)
    assert match is not None
    assert match.group(1) == version
