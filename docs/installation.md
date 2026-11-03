# Installation

`ggann` requires Python 3.12 or later.

ggann needs [`annplyr`](https://github.com/mdmanurung/annplyr)
`>=0.5,<0.6` and [`plotnine-extra`](https://github.com/mdmanurung/plotnine-extra)
`>=0.3.1,<0.4`. Neither release is on PyPI yet, so `pip install ggann` does not
resolve today. The package metadata does not use Git commit URLs, so install the
two dependencies from GitHub first, at the commits CI tests against:

```bash
pip install "annplyr @ git+https://github.com/mdmanurung/annplyr@ae2db5402c8b1c714bb12906ff82a40fbb0d8de2"
pip install "plotnine-extra @ git+https://github.com/mdmanurung/plotnine-extra@fcbd88f64760a2dd6e8200028a1a00909aa2c9f8"
pip install "ggann @ git+https://github.com/mdmanurung/ggann"
```

Once ggann 0.1.0 is released and both dependencies are on PyPI, one command
will do:

```bash
pip install ggann
```

The core install also pulls in [plotnine](https://plotnine.org/), Scanpy and
AnnData. All runtime requirements are declared as bounded package versions.

## Optional extras

| Extra | Functions | Backend |
|---|---|---|
| `density` | {func}`~ggann.plot_density` | [pyNebulosa](https://github.com/mdmanurung/pyNebulosa) |
| `upset` | {func}`~ggann.plot_upset` | [marsilea](https://marsilea.readthedocs.io/) |
| `heatmap` | {func}`~ggann.plot_clustermap`, {func}`~ggann.tidy_heatmap` | [PyComplexHeatmap](https://github.com/DingWB/PyComplexHeatmap) |
| `pseudobulk` | {func}`~ggann.pseudobulk` | [decoupler](https://decoupler-py.readthedocs.io/) |

Install one or more extras from the same distribution:

```bash
pip install "ggann[density,upset,heatmap,pseudobulk]"
```

Calling one of these functions without its backend raises an `ImportError` that
names the required extra.

## Install a source checkout

For release verification or development, install `annplyr` and `plotnine-extra`
as above, then:

```bash
git clone https://github.com/mdmanurung/ggann.git
cd ggann
python -m pip install -e ".[test,docs]"
```

Do not add sibling repositories to `PYTHONPATH` when validating a wheel. A clean
wheel-install test must resolve every runtime dependency from declared metadata.
