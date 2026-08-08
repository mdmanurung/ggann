# Installation

`ggann` requires Python 3.12 or later.

```bash
pip install git+https://github.com/mdmanurung/ggann
```

The core install includes [plotnine](https://plotnine.org/),
[`plotnine-extra`](https://github.com/mdmanurung/plotnine-extra), and
[`annplyr`](https://github.com/mdmanurung/annplyr).

## Optional extras

| Extra | Functions | Backend |
|---|---|---|
| `density` | {func}`~ggann.plot_density` | [pyNebulosa](https://github.com/mdmanurung/pyNebulosa) |
| `upset` | {func}`~ggann.plot_upset` | [marsilea](https://marsilea.readthedocs.io/) |
| `heatmap` | {func}`~ggann.plot_clustermap` | [PyComplexHeatmap](https://github.com/DingWB/PyComplexHeatmap) |
| `pseudobulk` | {func}`~ggann.pseudobulk` | [decoupler](https://decoupler-py.readthedocs.io/) |

Install one or more extras with the repository URL:

```bash
pip install "ggann[density,upset,heatmap,pseudobulk] @ git+https://github.com/mdmanurung/ggann"
```

Calling one of these functions without its backend raises an `ImportError` that
names the required extra.
