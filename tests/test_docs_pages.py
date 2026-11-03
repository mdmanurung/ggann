"""Every docs page shows at least one plot, and every plot it names exists."""

from __future__ import annotations

import re
from pathlib import Path

import pytest

DOCS = Path(__file__).resolve().parents[1] / "docs"
PAGES = ["index", "from_scanpy", "grammar", "publication", "heatmaps", "notes", "api"]
IMAGE = re.compile(r"(?:\]\(|\{image\}\s*|\{figure\}\s*|image::\s*)/?(images/[^)\s]+\.png)")


@pytest.mark.parametrize("page", PAGES)
def test_page_shows_existing_plots(page):
    text = (DOCS / f"{page}.md").read_text()
    images = IMAGE.findall(text)
    if page == "api":
        # api.md gets its plots from per-function images appended at build time
        images = [p.relative_to(DOCS).as_posix() for p in (DOCS / "images/api").glob("*.png")]
    assert images, f"{page}.md shows no plot"
    missing = [name for name in images if not (DOCS / name).is_file()]
    assert not missing, f"{page}.md names missing images: {missing}"


def test_site_has_only_the_listed_pages():
    pages = {p.stem for p in DOCS.glob("*.md")}
    assert pages == set(PAGES)
    assert not (DOCS / "vignettes").exists()
