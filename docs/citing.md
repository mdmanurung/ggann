# Citing ggann

If ggann made a figure in your paper, cite ggann and the tools that computed
what the figure shows. ggann draws data; the statistics and embeddings come
from the packages below.

## ggann

Cite the version you used. `CITATION.cff` in the repository holds the
metadata, and GitHub's "Cite this repository" button renders it as BibTeX.

```text
Manurung, M. ggann: a grammar of graphics for AnnData (version 0.1.0).
https://github.com/mdmanurung/ggann
```

Record the version with `importlib.metadata.version("ggann")` in your
methods, along with the versions of the packages below.

## Packages ggann always uses

| Package | What it does in ggann | Reference |
|---|---|---|
| AnnData | Holds the data every helper reads | Virshup I, Rybakov S, Theis FJ, Angerer P, Wolf FA. anndata: access and store annotated data matrices. *JOSS* 9, 4371 (2024). doi:10.21105/joss.04371 |
| Scanpy | Computes embeddings, clusters, `rank_genes_groups` and dendrograms that ggann plots | Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis. *Genome Biol* 19, 15 (2018). doi:10.1186/s13059-017-1382-0 |
| plotnine | Renders every ggplot-style figure | Kibirige H et al. plotnine: a grammar of graphics for Python. Zenodo, doi:10.5281/zenodo.1325308 (cite the version DOI you used) |
| SciPy | Runs the tests behind `stat_cor`, `stat_compare_means` and `stat_regline_equation` | Virtanen P et al. SciPy 1.0: fundamental algorithms for scientific computing in Python. *Nat Methods* 17, 261–272 (2020). doi:10.1038/s41592-019-0686-2 |

## Packages behind optional helpers

Cite these only when you use the helper.

| Helper | Package | Reference |
|---|---|---|
| `pseudobulk` | decoupler | Badia-i-Mompel P et al. decoupleR: ensemble of computational methods to infer biological activities from omics data. *Bioinformatics Advances* 2, vbac016 (2022). doi:10.1093/bioadv/vbac016 |
| `plot_clustermap` | PyComplexHeatmap | Ding W, Goldberg D, Zhou W. PyComplexHeatmap: a Python package to visualize multimodal genomics data. *iMeta* 2, e115 (2023). doi:10.1002/imt2.115 |
| `plot_density` | Nebulosa (via pynebulosa) | Alquicira-Hernandez J, Powell JE. Nebulosa recovers single-cell gene expression signals by kernel density estimation. *Bioinformatics* 37, 2485–2487 (2021). doi:10.1093/bioinformatics/btab003 |
| `plot_upset` | Marsilea | Zheng Y, Zheng Z, Rendeiro AF, Cheung E. Marsilea: an intuitive generalized paradigm for composable visualizations. *Genome Biol* 26, 5 (2025). doi:10.1186/s13059-024-03469-3 |

## What to report in methods

A figure is reproducible when the methods state:

- the expression source (`.X`, `.raw` or a named layer) and its normalisation;
- the statistic each panel shows (mean, fraction expressing, z-score);
- the test, its sidedness and the multiple-testing correction, for any
  p-value on a figure;
- for pseudobulk statistics, that samples or donors are the replicates;
- the `random_state` of any downsampling or jitter;
- how many cells were dropped for a missing group (ggann warns with the count).
