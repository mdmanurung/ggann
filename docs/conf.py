from __future__ import annotations

import os
import sys
from datetime import datetime
from importlib.metadata import PackageNotFoundError, metadata
from pathlib import Path

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE / "extensions"))

try:
    info = metadata("ggann")
    project = info["Name"]
    author = info["Author"] or "ggann developers"
    version = info["Version"]
    urls = dict(pu.split(", ") for pu in info.get_all("Project-URL") or [])
    repository_url = urls.get("Homepage", "https://github.com/mdmanurung/ggann")
except PackageNotFoundError:
    project = "ggann"
    author = "ggann developers"
    version = "0.1.0"
    repository_url = "https://github.com/mdmanurung/ggann"

copyright = f"{datetime.now():%Y}, {author}."
release = version

templates_path = ["_templates"]
nitpicky = False
needs_sphinx = "4.0"

html_context = {
    "display_github": True,
    "github_user": "mdmanurung",
    "github_repo": "ggann",
    "github_version": "main",
    "conf_py_path": "/docs/",
}

extensions = [
    "myst_nb",
    "sphinx_copybutton",
    "sphinx.ext.autodoc",
    "sphinx.ext.intersphinx",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx_autodoc_typehints",
    "sphinx_design",
    "IPython.sphinxext.ipython_console_highlighting",
    "sphinxext.opengraph",
    *[path.stem for path in (HERE / "extensions").glob("*.py")],
]

autosummary_generate = True
autodoc_member_order = "groupwise"
default_role = "literal"
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_use_rtype = True
myst_heading_anchors = 6
myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "dollarmath",
    "html_image",
    "html_admonition",
]
myst_url_schemes = ("http", "https", "mailto")
# Notebook execution remains off; ``executable_vignettes`` runs the six Python
# workflows in isolated processes before every build. They read the PBMC3K
# dataset from the ``data/`` cache that ``scripts/fetch_datasets.py`` populates,
# so the build itself reaches no network.
nb_execution_mode = "off"
typehints_defaults = "braces"
always_use_bars_union = True

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "myst-nb",
    ".ipynb": "myst-nb",
}

if os.environ.get("GGANN_DOCS_OFFLINE") == "1":
    intersphinx_mapping = {}
else:
    intersphinx_mapping = {
        "python": ("https://docs.python.org/3", None),
        "anndata": ("https://anndata.scverse.org/en/stable/", None),
        "scanpy": ("https://scanpy.readthedocs.io/en/stable/", None),
        "numpy": ("https://numpy.org/doc/stable/", None),
        "pandas": ("https://pandas.pydata.org/docs/", None),
        "plotnine": ("https://plotnine.org/", None),
    }

exclude_patterns = [
    "_build",
    # Older builds wrote every public re-export here, including third-party
    # objects whose upstream docstrings do not parse cleanly in Sphinx. Native
    # API pages now live under generated/native/.
    "generated/anngg.*",
    "generated/ggann.*",
    "Thumbs.db",
    ".DS_Store",
    "**.ipynb_checkpoints",
]

html_theme = "sphinx_book_theme"
html_static_path = ["_static"]
html_css_files = ["css/custom.css"]
html_title = project
html_theme_options = {
    "repository_url": repository_url,
    "use_repository_button": True,
    "path_to_docs": "docs/",
    "navigation_with_keys": False,
}

pygments_style = "default"
