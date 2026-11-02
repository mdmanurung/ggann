 12:25:48 up 72 days, 22:28,  0 users,  load average: 10.43, 17.17, 22.75
[1/3] extended/csr/x/embedding_categorical
[2/3] extended/csr/x/dotplot
[3/3] extended/csr/x/matrixplot
Wrote 3 matched cases to /tmp/claude-149488/-exports-para-lipg-hpc-mdmanurung--hobby-anngg/2fb8a426-3b10-48f4-a13c-58176b0988d3/scratchpad/perf/frozen.json
[1/3] extended/csr/x/embedding_categorical
[2/3] extended/csr/x/dotplot
[3/3] extended/csr/x/matrixplot
Wrote 3 matched cases to /tmp/claude-149488/-exports-para-lipg-hpc-mdmanurung--hobby-anngg/2fb8a426-3b10-48f4-a13c-58176b0988d3/scratchpad/perf/candidate.json
# Frozen ggann candidate regression check

Verdict: **FAIL**

Maximum permitted regression: 5.0%.
Baseline: `/tmp/claude-149488/-exports-para-lipg-hpc-mdmanurung--hobby-anngg/2fb8a426-3b10-48f4-a13c-58176b0988d3/scratchpad/perf/frozen.json` (`523a6c4ea03647e6cba6b652ce570c412f3dfcc21aea4b2f936e9bbddc23ad3e`).
Candidate: `/tmp/claude-149488/-exports-para-lipg-hpc-mdmanurung--hobby-anngg/2fb8a426-3b10-48f4-a13c-58176b0988d3/scratchpad/perf/candidate.json` (`cea6c19d0b894fa80fb6bf38cead504e4b5384610b99264f01d46cc29b0a49c5`).

| Workload | Metric | Frozen | Candidate | Candidate / frozen | Result |
|---|---|---:|---:|---:|:---:|
| dotplot | `preparation_median_seconds` | 38.805 ms | 38.969 ms | 1.004x | PASS |
| dotplot | `construction_median_seconds` | 43.094 ms | 44.088 ms | 1.023x | PASS |
| dotplot | `render_median_seconds` | 453.263 ms | 466.920 ms | 1.030x | PASS |
| dotplot | `end_to_end_median_seconds` | 686.129 ms | 700.983 ms | 1.022x | PASS |
| dotplot | `end_to_end_peak_rss_delta_bytes` | 13.383 MiB | 13.617 MiB | 1.018x | PASS |
| dotplot | `end_to_end_retained_after_gc_bytes` | 13.383 MiB | 13.766 MiB | 1.029x | PASS |
| embedding_categorical | `preparation_median_seconds` | 11.445 ms | 12.193 ms | 1.065x | FAIL |
| embedding_categorical | `construction_median_seconds` | 13.056 ms | 14.915 ms | 1.142x | FAIL |
| embedding_categorical | `render_median_seconds` | 354.354 ms | 373.337 ms | 1.054x | FAIL |
| embedding_categorical | `end_to_end_median_seconds` | 683.313 ms | 732.628 ms | 1.072x | FAIL |
| embedding_categorical | `end_to_end_peak_rss_delta_bytes` | 26.234 MiB | 26.820 MiB | 1.022x | PASS |
| embedding_categorical | `end_to_end_retained_after_gc_bytes` | 25.508 MiB | 25.184 MiB | 0.987x | PASS |
| matrixplot | `preparation_median_seconds` | 39.854 ms | 41.392 ms | 1.039x | PASS |
| matrixplot | `construction_median_seconds` | 43.002 ms | 45.655 ms | 1.062x | FAIL |
| matrixplot | `render_median_seconds` | 388.011 ms | 401.569 ms | 1.035x | PASS |
| matrixplot | `end_to_end_median_seconds` | 597.267 ms | 638.732 ms | 1.069x | FAIL |
| matrixplot | `end_to_end_peak_rss_delta_bytes` | 14.250 MiB | 14.324 MiB | 1.005x | PASS |
| matrixplot | `end_to_end_retained_after_gc_bytes` | 14.250 MiB | 14.324 MiB | 1.005x | PASS |

## Source trees

- Frozen: `fde527642c875ed653afda1bfca4cc049a3648db001e747b499969c8cc676f97`
- Candidate: `b100a7480f6dd5f901b8f5a7ab91c4a8d47d27933f86814b6731f05a0b7d21cb`
 13:24:10 up 72 days, 23:26,  0 users,  load average: 12.78, 12.75, 12.88
