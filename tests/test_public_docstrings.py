"""Release contract for first-party public API documentation."""

from __future__ import annotations

import inspect
import re

import ggann

_CALLABLE_SECTIONS = ("Parameters", "Returns", "Raises", "Notes", "Examples")


def _section(doc: str, name: str) -> str:
    match = re.search(
        rf"(?ms)^{name}\n-+\n(.*?)(?=^[A-Z][A-Za-z ]+\n-+\n|\Z)",
        doc,
    )
    assert match is not None, f"missing NumPy-style {name!r} section"
    return match.group(1)


def _first_party_public() -> dict[str, object]:
    return {
        name: getattr(ggann, name)
        for name in ggann.__all__
        if getattr(getattr(ggann, name), "__module__", "").startswith("ggann")
    }


def test_public_callables_have_complete_numpy_docstrings():
    for name, obj in _first_party_public().items():
        if not callable(obj):
            continue
        doc = inspect.getdoc(obj) or ""
        for heading in _CALLABLE_SECTIONS:
            _section(doc, heading)

        parameters = _section(doc, "Parameters")
        documented = {
            token.lstrip("*")
            for header in re.findall(r"(?m)^([^\n]+?)(?=\s*:|$)", parameters)
            for token in header.split(",")
            if token.strip()
            for token in [token.strip()]
        }
        expected = set(inspect.signature(obj).parameters)
        assert expected <= documented, (
            f"{name} is missing parameter docs for {sorted(expected - documented)}"
        )


def test_public_data_and_backend_boundaries_are_documented():
    sizes_doc = inspect.getdoc(ggann.sizes) or ""
    for heading in ("Attributes", "Notes", "Examples"):
        _section(sizes_doc, heading)

    assert "not a plotnine object" in inspect.getdoc(ggann.plot_clustermap).lower()
    assert "not a plotnine object" in inspect.getdoc(ggann.plot_upset).lower()
    assert ggann.scale_colour_expression is ggann.scale_color_expression
    assert ggann.scale_colour_celltype is ggann.scale_color_celltype
    assert ggann.scale_colour_obs is ggann.scale_color_obs
