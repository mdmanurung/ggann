"""ggann: a ggplot2-style plotting layer for scanpy AnnData objects.

The core grammar *is* plotnine. ``gganndata(adata) + aes(...) + geom_*()`` returns
a real :class:`plotnine.ggplot`, with aesthetics resolved from ``AnnData`` through
annplyr. High-level ``plot_*`` helpers reproduce scanpy's core figures, and
``plot_clustermap`` is a separate escape hatch for clustered heatmaps via
PyComplexHeatmap.
"""

from __future__ import annotations

from importlib.metadata import PackageNotFoundError, version

from plotnine_extra import (
    Beside,
    Stack,
    Wrap,
    geom_label_repel,
    geom_signif,
    geom_text_repel,
    ggsave,
    plot_annotation,
    plot_layout,
    stat_anova_test,
    stat_central_tendency,
    stat_compare_means,
    stat_cor,
    stat_kruskal_test,
    stat_pvalue_manual,
    stat_pwc,
    stat_regline_equation,
)

from ._matplotlib_backend import MatplotlibGGPlot
from ._palette import obs_colors, scale_color_obs, scale_colour_obs, scale_fill_obs
from ._resolve import embedding_coords, gene, obs, obsm
from .clustermap import plot_clustermap
from .composition import plot_proportions
from .correlation import plot_correlation
from .de import (
    plot_ma,
    plot_rank_genes_dotplot,
    plot_rank_genes_matrixplot,
    plot_volcano,
    rank_genes_df,
)
from .dendrogram import plot_dendrogram
from .density import plot_density
from .distributions import (
    plot_box,
    plot_expression_bar,
    plot_expression_line,
    plot_ridge,
    plot_sina,
    plot_stacked_violin,
    plot_violin,
)
from .grammar import aes, gganndata
from .layout import compose, tag_panels
from .markers import (
    plot_dotplot_grouped,
    plot_matrixplot_grouped,
    plot_tracksplot,
)
from .plots import (
    plot_dotplot,
    plot_embedding,
    plot_embedding_density,
    plot_features,
    plot_heatmap,
    plot_matrixplot,
)
from .pseudobulk import pseudobulk
from .qc import (
    plot_highest_expr_genes,
    plot_qc_scatter,
    plot_qc_violin,
    plot_variance_ratio,
)
from .theme import (
    reset_theme,
    scale_color_celltype,
    scale_color_expression,
    scale_colour_celltype,
    scale_colour_expression,
    scale_fill_celltype,
    scale_fill_expression,
    set_theme,
    sizes,
    theme_ggann,
)
from .upset import plot_upset

try:
    __version__ = version("ggann")
except PackageNotFoundError:  # pragma: no cover - not installed
    __version__ = "0.1.0"

__all__ = [
    "gganndata",
    "aes",
    "gene",
    "obs",
    "obsm",
    "embedding_coords",
    "MatplotlibGGPlot",
    "plot_embedding",
    "plot_embedding_density",
    "plot_features",
    "plot_dotplot",
    "plot_matrixplot",
    "plot_heatmap",
    "plot_violin",
    "plot_sina",
    "plot_ridge",
    "plot_clustermap",
    "plot_dendrogram",
    "plot_proportions",
    "plot_variance_ratio",
    "rank_genes_df",
    "plot_rank_genes_dotplot",
    "plot_rank_genes_matrixplot",
    "plot_volcano",
    "plot_ma",
    "plot_stacked_violin",
    "plot_tracksplot",
    "plot_dotplot_grouped",
    "plot_matrixplot_grouped",
    "plot_qc_violin",
    "plot_qc_scatter",
    "plot_highest_expr_genes",
    "plot_density",
    "plot_box",
    "plot_expression_bar",
    "plot_expression_line",
    "plot_correlation",
    "plot_upset",
    "theme_ggann",
    "set_theme",
    "reset_theme",
    "sizes",
    "scale_color_expression",
    "scale_colour_expression",
    "scale_fill_expression",
    "scale_color_celltype",
    "scale_colour_celltype",
    "scale_fill_celltype",
    "obs_colors",
    "scale_color_obs",
    "scale_colour_obs",
    "scale_fill_obs",
    "compose",
    "tag_panels",
    "Beside",
    "Stack",
    "Wrap",
    "plot_layout",
    "plot_annotation",
    "geom_text_repel",
    "geom_label_repel",
    "pseudobulk",
    "stat_compare_means",
    "stat_pwc",
    "stat_pvalue_manual",
    "stat_cor",
    "stat_regline_equation",
    "stat_anova_test",
    "stat_kruskal_test",
    "stat_central_tendency",
    "geom_signif",
    "ggsave",
    "__version__",
]
