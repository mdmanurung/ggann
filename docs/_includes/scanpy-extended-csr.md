Measured on 20,000 observations, 10,000 variables, 32 requested genes, and 32 groups (7 warm repetitions, CSR `.X`).

| Workload | ggann prep | Scanpy prep | ggann construct | Scanpy construct | ggann render/save | Scanpy render/save | ggann end to end | Scanpy end to end | E2E speedup |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Categorical embedding | 13.3 ms | 4.3 ms | 14.6 ms | 94.4 ms | 382.6 ms | 313.4 ms | 770.4 ms | 426.7 ms | 0.55x |
| Dotplot | 49.5 ms | 17.2 ms | 53.5 ms | 22.9 ms | 457.8 ms | 213.0 ms | 738.0 ms | 367.9 ms | 0.50x |
| Matrixplot | 51.6 ms | 12.0 ms | 55.7 ms | 15.2 ms | 413.9 ms | 217.1 ms | 656.7 ms | 346.5 ms | 0.53x |

Fresh-child memory probes exclude imports and fixture creation from the baseline:

| Workload | ggann peak RSS | Scanpy peak RSS | ggann retained RSS | Scanpy retained RSS | Peak ratio |
|---|---:|---:|---:|---:|---:|
| Categorical embedding | 29.09 MiB | 5.60 MiB | 28.17 MiB | 5.60 MiB | 5.19x |
| Dotplot | 17.00 MiB | 12.12 MiB | 17.00 MiB | 12.12 MiB | 1.40x |
| Matrixplot | 17.16 MiB | 7.73 MiB | 17.16 MiB | 7.73 MiB | 2.22x |

Speedup is Scanpy time divided by ggann time; values above 1 favour ggann. The recorded gate status is **FAIL**. All prepared payloads passed the benchmark's numeric equivalence checks and all input fingerprints were unchanged.

Construction is diagnostic because Scanpy eagerly creates different amounts of Matplotlib state. `render/save` measures PNG save from a materialized figure; end-to-end is the authoritative plotting comparison.

Recorded ggann source-tree SHA-256: `7e059e89a18746625fb9e3727e2b1f1d7d1ae22050d2b5dfa711832253bf4e0d`.
