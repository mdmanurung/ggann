# Installation

`ggann` requires Python 3.12 or later.

```bash
pip install ggann
```

The core install includes [plotnine](https://plotnine.org/),
[`plotnine-extra`](https://github.com/mdmanurung/plotnine-extra), and
[`annplyr`](https://github.com/mdmanurung/annplyr). A release is publishable only
after every required dependency is available from the configured package index;
the distribution metadata must not rely on Git commit URLs.

## Optional extras

| Extra | Functions | Backend |
|---|---|---|
| `density` | {func}`~ggann.plot_density` | [pyNebulosa](https://github.com/mdmanurung/pyNebulosa) |
| `upset` | {func}`~ggann.plot_upset` | [marsilea](https://marsilea.readthedocs.io/) |
| `heatmap` | {func}`~ggann.plot_clustermap` | [PyComplexHeatmap](https://github.com/DingWB/PyComplexHeatmap) |
| `pseudobulk` | {func}`~ggann.pseudobulk` | [decoupler](https://decoupler-py.readthedocs.io/) |

Install one or more extras from the same distribution:

```bash
pip install "ggann[density,upset,heatmap,pseudobulk]"
```

Calling one of these functions without its backend raises an `ImportError` that
names the required extra.

## Install a source checkout

For release verification or development:

```bash
git clone https://github.com/mdmanurung/ggann.git
cd ggann
python -m pip install -e ".[test,docs]"
```

Do not add sibling repositories to `PYTHONPATH` when validating a wheel. A clean
wheel-install test must resolve every runtime dependency from declared metadata.
