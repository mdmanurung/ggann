# Frozen ggann candidate regression check

Verdict: **PASS**

Maximum permitted regression: 5.0%.
Baseline: `publication-baseline-30e3359.json` (`a5f5187263e8ea448f39bd160c77047928f2724013a0cca59998ac76f93d6fd9`).
Candidate: `publication-candidate.json` (`2ec917890297c340dff862f1a3077e873b3b11ee0c3ccd41a127939ced868587`).

| Workload | Metric | Frozen | Candidate | Candidate / frozen | Result |
|---|---|---:|---:|---:|:---:|
| dotplot | `preparation_median_seconds` | 39.830 ms | 39.945 ms | 1.003x | PASS |
| dotplot | `construction_median_seconds` | 43.781 ms | 43.854 ms | 1.002x | PASS |
| dotplot | `render_median_seconds` | 516.347 ms | 516.490 ms | 1.000x | PASS |
| dotplot | `end_to_end_median_seconds` | 770.675 ms | 767.436 ms | 0.996x | PASS |
| dotplot | `end_to_end_peak_rss_delta_bytes` | 14.133 MiB | 13.168 MiB | 0.932x | PASS |
| dotplot | `end_to_end_retained_after_gc_bytes` | 14.133 MiB | 13.168 MiB | 0.932x | PASS |
| embedding_categorical | `preparation_median_seconds` | 13.794 ms | 13.279 ms | 0.963x | PASS |
| embedding_categorical | `construction_median_seconds` | 15.053 ms | 14.357 ms | 0.954x | PASS |
| embedding_categorical | `render_median_seconds` | 400.124 ms | 399.288 ms | 0.998x | PASS |
| embedding_categorical | `end_to_end_median_seconds` | 769.493 ms | 764.310 ms | 0.993x | PASS |
| embedding_categorical | `end_to_end_peak_rss_delta_bytes` | 25.461 MiB | 25.512 MiB | 1.002x | PASS |
| embedding_categorical | `end_to_end_retained_after_gc_bytes` | 24.582 MiB | 23.598 MiB | 0.960x | PASS |
| matrixplot | `preparation_median_seconds` | 40.154 ms | 39.858 ms | 0.993x | PASS |
| matrixplot | `construction_median_seconds` | 44.550 ms | 44.969 ms | 1.009x | PASS |
| matrixplot | `render_median_seconds` | 464.766 ms | 461.646 ms | 0.993x | PASS |
| matrixplot | `end_to_end_median_seconds` | 698.835 ms | 699.023 ms | 1.000x | PASS |
| matrixplot | `end_to_end_peak_rss_delta_bytes` | 14.492 MiB | 13.531 MiB | 0.934x | PASS |
| matrixplot | `end_to_end_retained_after_gc_bytes` | 14.492 MiB | 13.531 MiB | 0.934x | PASS |

## Source trees

- Frozen: `fde527642c875ed653afda1bfca4cc049a3648db001e747b499969c8cc676f97`
- Candidate: `f834d9443a244a6d9472280b6df813789693b1467ca15510f22046482b7a7ab2`
