"""Append a rendered example image to each API-reference page that has one.

Images are pre-rendered by ``examples/api_examples.py`` into
``docs/images/api/<fullname>.png`` and committed. This hook adds an ``Examples``
section pointing at the image for any documented object that has one; objects
without an image (scales, grammar accessors, ``rank_genes_df``) are untouched, so
there are no broken images and the source docstrings stay clean.
"""

from __future__ import annotations

import os

from sphinx.application import Sphinx


def _process_docstring(app, what, name, obj, options, lines):
    if what not in ("function", "method"):
        return
    img = os.path.join(app.srcdir, "images", "api", f"{name}.png")
    if not os.path.exists(img):
        return
    lines += [
        "",
        ".. rubric:: Rendered output",
        "",
        f".. image:: /images/api/{name}.png",
        "   :alt: Example output of " + name,
        "   :width: 90%",
        "",
    ]


def setup(app: Sphinx) -> dict:
    app.connect("autodoc-process-docstring", _process_docstring)
    return {"parallel_read_safe": True, "parallel_write_safe": True}
