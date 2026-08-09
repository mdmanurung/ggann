# Downstream artifact invariance

| Artifact | Result | Detail |
|---|:---:|---|
| provenance and parameters | pass | exact comparison; source digests recorded separately |
| `resolve` data | pass | equal within rtol 1e-6, atol 1e-7 |
| `aggregate` data | pass | equal within rtol 1e-6, atol 1e-7 |
| `plot_embedding` data | pass | equal within rtol 1e-6, atol 1e-7 |
| `plot_dotplot` data | pass | equal within rtol 1e-6, atol 1e-7 |
| `plot_matrixplot` data | pass | equal within rtol 1e-6, atol 1e-7 |
| `plot_heatmap` data | pass | equal within rtol 1e-6, atol 1e-7 |
| plot mappings | pass | exact comparison |
| plot labels and scales | pass | exact comparison |
| `embedding` artists | pass | 300 points; 0 tiles; 0 rasterized collections; guide and axis text exact |
| `dotplot` artists | pass | 40 points; 0 tiles; 0 rasterized collections; guide and axis text exact |
| `matrixplot` artists | pass | 0 points; 40 tiles; 0 rasterized collections; guide and axis text exact |
| `heatmap` artists | pass | 0 points; 1244 tiles; 0 rasterized collections; guide and axis text exact |
| `embedding` SVG | pass | vector collections; no embedded image |
| `dotplot` SVG | pass | vector collections; no embedded image |
| `matrixplot` SVG | pass | vector collections; no embedded image |
| `embedding` image | pass | max 0; mean 0; PSNR inf dB |
| `dotplot` image | pass | max 0; mean 0; PSNR inf dB |
| `matrixplot` image | pass | max 0; mean 0; PSNR inf dB |
| `heatmap` image | pass | max 0; mean 0; PSNR inf dB |
