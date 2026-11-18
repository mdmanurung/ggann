# Shapes

`ggann.shapes` adds layers for shapes plotnine does not draw. This page covers curves. Every layer is a plain plotnine layer, so you add it with `+` and combine it with any geom, scale, facet or theme.

All examples come from `examples/docs_figures/shapes.py`.

## Curves

`geom_curve` draws one bent line for each row of a data frame with `x`, `y`, `xend` and `yend`. `stat_curve` does the work: it replaces each row by `n` points on a quadratic Bezier curve and adds a `t` column from 0 to 1.

```python
import pandas as pd
from plotnine import aes, arrow, coord_equal, ggplot
import ggann as ag

df = pd.DataFrame({"x": [0.0], "y": [0.0], "xend": [4.0], "yend": [0.0]})
base = ggplot(df, aes("x", "y", xend="xend", yend="yend")) + coord_equal()
plot = base + ag.geom_curve(curvature=0.5, arrow=arrow(length=0.15))
```

The figure draws one segment from `(0, 0)` to `(4, 0)` three times, with `curvature` -0.5, 0 and 0.5. Each curve travels left to right.

![One segment drawn as three arrows: a blue curve bending up for curvature -0.5, a straight grey line for 0 and a red curve bending down for 0.5](images/shapes/curve.png)

Positive `curvature` bends right of the direction of travel, negative bends left, and 0 is a straight line. The bulge at the middle is `|curvature| / 2` of the segment length. The defaults are `curvature=0.5` and `n=50`. `n` must be an integer of at least 2, and `curvature` must be finite. A bad value raises `ValueError` when the plot is built.

Colour, size, alpha and linetype map per row and stay constant along a curve. `arrow`, `lineend` and `linejoin` go to the path geom.

### Limits

- **Curvature lives in data units.** The bend is computed in transformed scale units, so the shape is exact only when x and y share units, as under `coord_equal()`. A segment from `(0, 0)` to `(1, 1000)` at the default curvature reaches x = 250.5 and stretches the x scale about 250-fold. Under `scale_x_log10` the curve bends in log units.
- **Not for `coord_radial`.** There the path bends in (theta, r) space and can push r below 0. Curves under `coord_radial` are unsupported.
- **`group` is set by the stat.** It numbers the surviving rows 1 to k in each panel and replaces any mapped `group`.
- **Missing endpoints are dropped.** A row with a missing or non-finite `x`, `y`, `xend` or `yend` is removed with a warning. `na_rm=True` silences it.
- **Registry names are global.** plotnine keeps one registry keyed by class name. After `import ggann`, `geom_path(stat="curve")` resolves to `stat_curve` in every plot in the process. A same-named class from another package would replace ggann's silently, and a canary test checks that none exists.

There is no `ggcurve` constructor. Use the layers.
