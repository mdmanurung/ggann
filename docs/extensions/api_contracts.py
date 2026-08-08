"""Add concise usage contracts to first-party API pages."""

from __future__ import annotations

from sphinx.application import Sphinx
from sphinx.util import logging

LOGGER = logging.getLogger(__name__)

_PBMC_SETUP = """\
import scanpy as sc
import ggann as ag

adata = sc.datasets.pbmc68k_reduced()"""

_DE_SETUP = (
    _PBMC_SETUP
    + """
sc.tl.rank_genes_groups(adata, "bulk_labels", method="wilcoxon")"""
)

_EXPRESSION_ERRORS = (
    "Raises ``KeyError`` when a requested gene, observation column, or layer is "
    "missing. Incompatible expression sources and invalid ordering, scaling, or "
    "downsampling options raise ``ValueError``."
)
_EMBEDDING_ERRORS = (
    "A missing embedding or colour field raises ``KeyError``. An embedding with "
    "fewer than two coordinates, or an invalid downsampling value, raises ``ValueError``."
)
_PLOTNINE_ERRORS = (
    "There is no additional ggann validation; plotnine reports unsupported layer, "
    "scale, or theme arguments."
)


def _pbmc(call: str) -> str:
    return _PBMC_SETUP + "\n\n" + call


def _de(call: str) -> str:
    return _DE_SETUP + "\n\n" + call


def _contract(returns: str, errors: str, example: str) -> dict[str, str]:
    return {"returns": returns, "errors": errors, "example": example}


_CONTRACTS: dict[str, dict[str, str]] = {
    "ggann.gganndata": _contract(
        "A regular ``plotnine.ggplot`` built from the resolved per-cell data.",
        "An explicit missing gene, observation column, embedding, or layer raises "
        "``KeyError``. Combining ``layer=`` with ``use_raw=True`` raises ``ValueError``.",
        _pbmc(
            "from plotnine import geom_point\n\n"
            "plot = ag.gganndata(\n"
            '    adata, ag.aes("UMAP_1", "UMAP_2", color="bulk_labels")\n'
            ") + geom_point()"
        ),
    ),
    "ggann.gene": _contract(
        "A source-qualified gene reference for an aesthetic mapping.",
        "Construction does not inspect ``adata``. A missing gene or incompatible "
        "source is reported when the reference is resolved.",
        _pbmc(
            'mapping = ag.aes(color=ag.gene("CD3D", use_raw=True))\n'
            "plot = ag.gganndata(adata, mapping)"
        ),
    ),
    "ggann.obs": _contract(
        "A source-qualified observation-column reference.",
        "Construction does not inspect ``adata``. A missing observation column is "
        "reported when the reference is resolved.",
        _pbmc(
            'mapping = ag.aes(color=ag.obs("bulk_labels"))\n'
            "plot = ag.gganndata(adata, mapping)"
        ),
    ),
    "ggann.obsm": _contract(
        "A reference to one zero-based coordinate of an ``obsm`` matrix.",
        "A non-integer coordinate fails during construction. A missing basis or "
        "out-of-range coordinate is reported when the reference is resolved.",
        _pbmc(
            'mapping = ag.aes(x=ag.obsm("umap", 0), y=ag.obsm("umap", 1))\n'
            "plot = ag.gganndata(adata, mapping)"
        ),
    ),
    "ggann.embedding_coords": _contract(
        "A ``pandas.DataFrame`` indexed like ``adata.obs``.",
        "A missing basis raises ``KeyError``; a negative ``n`` raises ``ValueError``.",
        _pbmc('coordinates = ag.embedding_coords(adata, "umap", n=2)'),
    ),
}


_GGPLOT_CALLS = {
    "plot_embedding": 'plot = ag.plot_embedding(adata, "umap", color="bulk_labels")',
    "plot_embedding_density": (
        "plot = ag.plot_embedding_density(\n"
        '    adata, "umap", group_by="bulk_labels"\n'
        ")"
    ),
    "plot_features": 'plot = ag.plot_features(adata, ["CD3D", "NKG7"])',
    "plot_density": 'plot = ag.plot_density(adata, ["CD3D", "NKG7"])',
    "plot_dotplot": (
        "plot = ag.plot_dotplot(\n"
        '    adata, ["CD3D", "NKG7"], group_by="bulk_labels"\n'
        ")"
    ),
    "plot_dotplot_grouped": (
        'genes = {"T cells": ["CD3D"], "NK cells": ["NKG7"]}\n'
        'plot = ag.plot_dotplot_grouped(adata, genes, group_by="bulk_labels")'
    ),
    "plot_matrixplot": (
        "plot = ag.plot_matrixplot(\n"
        '    adata, ["CD3D", "NKG7"], group_by="bulk_labels"\n'
        ")"
    ),
    "plot_matrixplot_grouped": (
        'genes = {"T cells": ["CD3D"], "NK cells": ["NKG7"]}\n'
        'plot = ag.plot_matrixplot_grouped(adata, genes, group_by="bulk_labels")'
    ),
    "plot_heatmap": (
        "plot = ag.plot_heatmap(\n"
        '    adata, ["CD3D", "NKG7"], group_by="bulk_labels",\n'
        '    standard_scale="var",\n'
        ")"
    ),
    "plot_violin": (
        'plot = ag.plot_violin(\n    adata, ["CD3D", "NKG7"], group_by="bulk_labels"\n)'
    ),
    "plot_ridge": (
        'plot = ag.plot_ridge(\n    adata, ["CD3D", "NKG7"], group_by="bulk_labels"\n)'
    ),
    "plot_stacked_violin": (
        "plot = ag.plot_stacked_violin(\n"
        '    adata, ["CD3D", "NKG7"], group_by="bulk_labels"\n'
        ")"
    ),
    "plot_tracksplot": (
        "plot = ag.plot_tracksplot(\n"
        '    adata, ["CD3D", "NKG7"], group_by="bulk_labels"\n'
        ")"
    ),
    "plot_dendrogram": 'plot = ag.plot_dendrogram(adata, group_by="bulk_labels")',
    "plot_box": (
        'plot = ag.plot_box(\n    adata, ["CD3D", "NKG7"], group_by="bulk_labels"\n)'
    ),
    "plot_sina": (
        'plot = ag.plot_sina(\n    adata, ["CD3D", "NKG7"], group_by="bulk_labels"\n)'
    ),
    "plot_expression_bar": (
        "plot = ag.plot_expression_bar(\n"
        '    adata, ["CD3D", "NKG7"], group_by="bulk_labels"\n'
        ")"
    ),
    "plot_expression_line": (
        "plot = ag.plot_expression_line(\n"
        '    adata, ["CD3D"], x="phase", group_by="bulk_labels"\n'
        ")"
    ),
    "plot_proportions": (
        "plot = ag.plot_proportions(\n"
        '    adata, group_by="bulk_labels", split_by="phase"\n'
        ")"
    ),
    "plot_correlation": 'plot = ag.plot_correlation(adata, group_by="bulk_labels")',
    "plot_qc_violin": (
        'plot = ag.plot_qc_violin(\n    adata, metrics=["n_genes", "percent_mito"]\n)'
    ),
    "plot_qc_scatter": ('plot = ag.plot_qc_scatter(adata, x="n_counts", y="n_genes")'),
    "plot_highest_expr_genes": (
        "plot = ag.plot_highest_expr_genes(adata, n=10, use_raw=True)"
    ),
    "plot_variance_ratio": "plot = ag.plot_variance_ratio(adata, n_pcs=20)",
}

_GGPLOT_ERRORS = {name: _EXPRESSION_ERRORS for name in _GGPLOT_CALLS}
for _name in ("plot_embedding", "plot_embedding_density", "plot_features"):
    _GGPLOT_ERRORS[_name] = _EMBEDDING_ERRORS
_GGPLOT_ERRORS.update(
    {
        "plot_features": (
            "A missing embedding raises ``KeyError``. If no requested gene or numeric "
            "observation field resolves, or downsampling is invalid, the function "
            "raises ``ValueError``."
        ),
        "plot_density": (
            "Raises ``ImportError`` when the ``density`` extra is unavailable, "
            "``KeyError`` for missing features or embeddings, and ``TypeError`` "
            "for non-numeric features."
        ),
        "plot_dendrogram": (
            "An unsupported orientation raises ``ValueError``. Missing or invalid "
            "grouping data are reported by scanpy when the tree is computed."
        ),
        "plot_proportions": (
            "Missing grouping columns raise ``KeyError``. Unsupported ``kind`` values "
            "or trend/area plots without ``split_by`` raise ``ValueError``."
        ),
        "plot_correlation": (
            "Missing genes, grouping columns, or layers raise ``KeyError``. Fewer than "
            "two selected genes and unsupported correlation methods raise ``ValueError``."
        ),
        "plot_qc_violin": (
            "Missing requested metrics raise ``KeyError``. If no default QC metrics "
            "exist, the function raises ``ValueError``."
        ),
        "plot_qc_scatter": (
            "A field that cannot be resolved from observations or expression raises "
            "``KeyError``; incompatible expression sources raise ``ValueError``."
        ),
        "plot_highest_expr_genes": (
            "A non-positive ``n`` or incompatible expression sources raise "
            "``ValueError``; a missing layer raises ``KeyError``."
        ),
        "plot_variance_ratio": (
            "Missing PCA variance information raises ``KeyError``."
        ),
    }
)

for _name, _call in _GGPLOT_CALLS.items():
    _CONTRACTS[f"ggann.{_name}"] = _contract(
        "A composable ``plotnine.ggplot`` object.",
        _GGPLOT_ERRORS[_name],
        _pbmc(_call),
    )

for _name, _call in {
    "plot_rank_genes_dotplot": (
        "plot = ag.plot_rank_genes_dotplot(\n"
        '    adata, n_genes=3, group_by="bulk_labels"\n'
        ")"
    ),
    "plot_rank_genes_matrixplot": (
        "plot = ag.plot_rank_genes_matrixplot(\n"
        '    adata, n_genes=3, group_by="bulk_labels"\n'
        ")"
    ),
    "plot_volcano": 'plot = ag.plot_volcano(adata, group="CD56+ NK")',
}.items():
    _errors = (
        "A missing rank-genes result raises ``KeyError``. Results without p-values "
        "or fold changes raise ``ValueError`` when those fields are required."
    )
    _CONTRACTS[f"ggann.{_name}"] = _contract(
        "A composable ``plotnine.ggplot`` object.", _errors, _de(_call)
    )

_CONTRACTS.update(
    {
        "ggann.rank_genes_df": _contract(
            "A tidy ``pandas.DataFrame`` with one row per ranked gene and group.",
            "A missing ``adata.uns[key]`` rank-genes result raises ``KeyError``; "
            "scanpy validates group names and filters.",
            _de('table = ag.rank_genes_df(adata, group="CD56+ NK")'),
        ),
        "ggann.plot_ma": _contract(
            "A composable ``plotnine.ggplot`` object.",
            "Missing mean, fold-change, or adjusted-p-value columns raise ``KeyError``.",
            """\
import pandas as pd
import ggann as ag

results = pd.DataFrame(
    {
        "baseMean": [10.0, 40.0, 100.0],
        "log2FoldChange": [-1.2, 0.1, 1.5],
        "padj": [0.01, 0.8, 0.02],
    },
    index=["gene_a", "gene_b", "gene_c"],
)
plot = ag.plot_ma(results)""",
        ),
        "ggann.plot_clustermap": _contract(
            "A ``PyComplexHeatmap.ClusterMapPlotter`` instance.",
            "Raises ``ImportError`` when the ``heatmap`` extra is unavailable, "
            "``KeyError`` for missing data, and ``ValueError`` when "
            "``standard_scale`` and ``z_score`` are both set.",
            _pbmc(
                "clustered = ag.plot_clustermap(\n"
                '    adata, ["CD3D", "NKG7"], group_by="bulk_labels"\n'
                ")"
            ),
        ),
        "ggann.plot_upset": _contract(
            "A marsilea ``Upset`` object, rendered immediately by default.",
            "Raises ``ImportError`` when the ``upset`` extra is unavailable, "
            "``TypeError`` when ``sets`` is not a mapping, and ``ValueError`` for "
            "fewer than two sets.",
            """\
import ggann as ag

sets = {"T cells": {"CD3D", "IL7R"}, "NK cells": {"NKG7", "GNLY"}}
upset = ag.plot_upset(sets)""",
        ),
        "ggann.pseudobulk": _contract(
            "A new ``anndata.AnnData`` with one observation per retained sample/group profile.",
            "Raises ``ImportError`` when the ``pseudobulk`` extra is unavailable and "
            "``ValueError`` when ``layer=`` is combined with ``use_raw=True``. "
            "decoupler validates count data and observation columns.",
            """\
import anndata as ad
import ggann as ag
import numpy as np
import pandas as pd

counts = ad.AnnData(
    X=np.array([[1, 0], [2, 1], [0, 3], [1, 2]], dtype=int),
    obs=pd.DataFrame(
        {
            "donor": ["d1", "d1", "d2", "d2"],
            "cell_type": ["T", "T", "T", "T"],
        },
        index=["cell1", "cell2", "cell3", "cell4"],
    ),
)
pb = ag.pseudobulk(
    counts, sample_col="donor", group_by="cell_type", min_cells=1
)""",
        ),
    }
)


_CONTRACTS.update(
    {
        "ggann.theme_ggann": _contract(
            "A plotnine theme object.",
            _PLOTNINE_ERRORS,
            "import ggann as ag\n\ntheme = ag.theme_ggann(base_size=9)",
        ),
        "ggann.set_theme": _contract(
            "The configured plotnine theme object.",
            "Invalid font or theme arguments are reported by plotnine. With "
            "``register=True``, this function intentionally changes global plotnine state.",
            "import ggann as ag\n\ntheme = ag.set_theme(base_size=9, register=False)",
        ),
        "ggann.reset_theme": _contract(
            "``None``.",
            "There are no ggann-specific error conditions.",
            "import ggann as ag\n\nag.set_theme()\nag.reset_theme()",
        ),
        "ggann.sizes": _contract(
            "The shared ggann font-size scale, measured in points.",
            "Calling ``geom`` with a non-numeric value raises ``TypeError``.",
            "import ggann as ag\n\npoint_size = ag.sizes.geom(ag.sizes.small)",
        ),
        "ggann.obs_colors": _contract(
            "A category-to-colour ``dict``, or ``None`` when no stored palette is usable.",
            "A missing or non-categorical column returns ``None`` rather than raising.",
            _pbmc('palette = ag.obs_colors(adata, "bulk_labels")'),
        ),
        "ggann.scale_color_obs": _contract(
            "A discrete plotnine colour scale.",
            _PLOTNINE_ERRORS,
            _pbmc('scale = ag.scale_color_obs(adata, "bulk_labels")'),
        ),
        "ggann.scale_fill_obs": _contract(
            "A discrete plotnine fill scale.",
            _PLOTNINE_ERRORS,
            _pbmc('scale = ag.scale_fill_obs(adata, "bulk_labels")'),
        ),
        "ggann.scale_color_expression": _contract(
            "A continuous plotnine colour scale.",
            _PLOTNINE_ERRORS,
            'import ggann as ag\n\nscale = ag.scale_color_expression(cmap="magma")',
        ),
        "ggann.scale_fill_expression": _contract(
            "A continuous plotnine fill scale.",
            _PLOTNINE_ERRORS,
            'import ggann as ag\n\nscale = ag.scale_fill_expression(cmap="viridis")',
        ),
        "ggann.scale_color_celltype": _contract(
            "A discrete plotnine colour scale.",
            _PLOTNINE_ERRORS,
            "import ggann as ag\n\nscale = ag.scale_color_celltype()",
        ),
        "ggann.scale_fill_celltype": _contract(
            "A discrete plotnine fill scale.",
            _PLOTNINE_ERRORS,
            "import ggann as ag\n\nscale = ag.scale_fill_celltype()",
        ),
        "ggann.compose": _contract(
            "A plotnine composition object with a ``save`` method.",
            "An empty panel sequence or unsupported tag level raises ``ValueError``.",
            _pbmc(
                'left = ag.plot_embedding(adata, "umap", color="bulk_labels")\n'
                "right = ag.plot_dotplot(\n"
                '    adata, ["CD3D", "NKG7"], group_by="bulk_labels"\n'
                ")\n"
                "figure = ag.compose([left, right], ncol=2)"
            ),
        ),
        "ggann.tag_panels": _contract(
            "A new list of tagged ``plotnine.ggplot`` objects.",
            "An unsupported ``levels`` value raises ``ValueError``.",
            _pbmc(
                "plots = [\n"
                '    ag.plot_embedding(adata, "umap", color="bulk_labels"),\n'
                '    ag.plot_embedding(adata, "umap", color="CD3D"),\n'
                "]\n"
                'tagged = ag.tag_panels(plots, levels="A")'
            ),
        ),
    }
)


def _process_docstring(app, what, name, obj, options, lines):
    if not name.startswith("ggann.") or what not in {"function", "data"}:
        return
    contract = _CONTRACTS.get(name)
    if contract is None:
        LOGGER.warning(
            "First-party API page %s has no return/error/example contract.",
            name,
            type="ggann.api_contract",
        )
        return

    lines += [
        "",
        ".. rubric:: Return value",
        "",
        contract["returns"],
        "",
        ".. rubric:: Notable errors",
        "",
        contract["errors"],
        "",
        ".. rubric:: Minimal example",
        "",
        ".. code-block:: python",
        "",
    ]
    lines.extend(
        f"   {line}" if line else "   " for line in contract["example"].splitlines()
    )
    lines.append("")


def setup(app: Sphinx) -> dict:
    app.connect("autodoc-process-docstring", _process_docstring)
    return {"parallel_read_safe": True, "parallel_write_safe": True}
