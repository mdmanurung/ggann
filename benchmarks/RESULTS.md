# Benchmark results

These measurements were collected on 2026-08-08 from the same checkout and
Python environment. The baseline used a copy of `src/ggann` taken before the
changes. Both runs used Python 3.14.5, anndata 0.12.19, annplyr 0.2.0, NumPy
2.4.6, pandas 2.3.3, plotnine 0.15.7, and SciPy 1.16.3. BLAS and OpenMP thread
counts were pinned to one.

The runner authenticated the baseline source as SHA-256
`255b30c3ea266de3d2280f3329f7302d3e48299e83f7e4301a8b4fd50ffd7b4d`
and the candidate source as
`0aa369844b38b7192e0c3007d57d335df30435cc98bb36eb40b0d8c903caa6d4`.

The standard fixture contains 4,000 observations, 2,000 variables, 16 requested
genes, and 16 groups. Sparse fixtures have 0.5% density. Each repeated result is
the median of five calls. Peak RSS is the largest increase over the pre-call RSS
sampled at 1 ms intervals.

| Format and workload | Median before | Median after | Time change | Peak RSS before | Peak RSS after | RSS change |
|---|---:|---:|---:|---:|---:|---:|
| CSC aggregate | 3.945 s | 203.8 ms | -94.8% | 22.9 MiB | 3.5 MiB | -84.7% |
| CSC resolve | 598.7 ms | 19.4 ms | -96.8% | 13.6 MiB | 3.3 MiB | -75.6% |
| CSC tidy | 649.5 ms | 43.6 ms | -93.3% | 15.6 MiB | 7.5 MiB | -51.8% |
| CSR aggregate | 3.956 s | 200.1 ms | -94.9% | 23.0 MiB | 3.6 MiB | -84.6% |
| CSR highest-expressed-gene preparation | 1.688 s | 15.6 ms | -99.1% | 197.4 MiB | 5.3 MiB | -97.3% |
| CSR resolve | 589.5 ms | 19.6 ms | -96.7% | 13.4 MiB | 3.0 MiB | -77.3% |
| CSR tidy | 608.4 ms | 46.0 ms | -92.4% | 15.8 MiB | 7.9 MiB | -49.8% |
| Dense aggregate | 326.8 ms | 43.3 ms | -86.8% | 95.3 MiB | 4.4 MiB | -95.4% |
| Dense resolve | 34.9 ms | 12.3 ms | -64.8% | 32.6 MiB | 3.0 MiB | -90.8% |
| Dense tidy | 50.2 ms | 34.4 ms | -31.5% | 32.5 MiB | 7.6 MiB | -76.5% |

No repeated runtime or peak-RSS case exceeded the 5% regression threshold.

Profiling the original 4,000 x 2,000 CSR aggregation placed 7.17 of 7.20 seconds
inside two annplyr summaries. Frame construction accounted for 2.44 seconds and
the complex fraction-expression group-by accounted for 4.21 seconds. The new
path projects requested variables before annplyr and summarizes means over a
bounded boolean expression projection, so both reductions use annplyr's simple
group-by path.

The highest-expressed-gene workload must rank every variable before selecting
the top genes. Its sparse-aware matrix reduction avoids the baseline's full
cells-by-variables DataFrame and is the measured reason this whole-matrix path
uses the centralized direct-matrix accessor. Its prepared data has identical
gene order and agrees within `rtol=1e-6`, `atol=1e-7` (maximum absolute
difference 2.73e-6; mean absolute difference 3.27e-9).

Exact output fingerprints changed where the optimized path retains float32
values instead of widening them, or changes the order of
mathematically equivalent sparse reductions. A provenance-authenticated
tolerance comparison on `scanpy.datasets.pbmc68k_reduced` found equal resolved
data, group summaries, prepared data, aesthetic mappings, labels, scales, and
categorical order. Rendered embedding, dotplot, matrixplot, and heatmap images
were pixel-identical (maximum and mean channel difference both zero). The
matrixplot intentionally omits the unused fraction-expression column from its
prepared data.

## Expression sources and rendering

The same standard CSR fixture was run with three repeats for named layers,
`.raw`, mixed per-aesthetic sources, plot preparation, and rendering.

| Workload | Median before | Median after | Time change | Peak RSS before | Peak RSS after |
|---|---:|---:|---:|---:|---:|
| Layer resolution | 615.0 ms | 21.1 ms | -96.6% | 13.4 MiB | 3.1 MiB |
| Raw resolution | 632.5 ms | 17.4 ms | -97.3% | 13.4 MiB | 3.0 MiB |
| Mixed-source resolution | 1.786 s | 26.8 ms | -98.5% | 12.2 MiB | 2.8 MiB |
| Layer aggregation | 3.796 s | 211.8 ms | -94.4% | 23.0 MiB | 4.0 MiB |
| Raw aggregation | 3.991 s | 200.3 ms | -95.0% | 21.9 MiB | 3.7 MiB |
| Dotplot preparation | 3.772 s | 231.3 ms | -93.9% | 23.1 MiB | 4.7 MiB |
| Dotplot rendering | 4.284 s | 464.1 ms | -89.2% | 28.3 MiB | 16.6 MiB |
| Embedding preparation | 19.2 ms | 17.2 ms | -10.3% | 4.3 MiB | 3.1 MiB |
| Embedding rendering | 122.9 ms | 121.7 ms | -1.0% | 14.5 MiB | 13.8 MiB |

## Intermediate sizes

On the standard fixture, the selected 4,000 x 16 expression block occupied
256,000 bytes for dense input, 18,596 bytes for CSR, and 2,660 bytes for CSC.
The prepared resolve frame occupied 628,692 bytes, the inherently long tidy
frame 4,145,512 bytes, and the 32-group aggregate frame 64,774 bytes. Isolated
scaling sweeps were not repeated for this final digest, so historical scaling
measurements are not attributed to it.

The benchmark inputs are synthetic and the RSS sampler can miss allocations
shorter than 1 ms. The first table excludes rendering; preparation and rendering
are reported separately above.
