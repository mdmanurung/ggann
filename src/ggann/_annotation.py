"""Cell-aware, contrast-aware vector annotations for matrix plot families."""

from __future__ import annotations

from typing import Any

import pandas as pd
from matplotlib.colors import to_rgb
from matplotlib.text import Text
from plotnine.geoms.geom_text import geom_text


def _relative_luminance(colour: object) -> float:
    try:
        channels = to_rgb(colour)
    except (TypeError, ValueError):
        return 1.0

    def linear(channel: float) -> float:
        return channel / 12.92 if channel <= 0.04045 else ((channel + 0.055) / 1.055) ** 2.4

    red, green, blue = (linear(channel) for channel in channels)
    return 0.2126 * red + 0.7152 * green + 0.0722 * blue


class _CellText(Text):
    """Text that decides visibility from the final rendered panel dimensions."""

    def __init__(self, *args, min_cell_pt: float, n_x: int, n_y: int, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self._min_cell_pt = min_cell_pt
        self._n_x = max(1, n_x)
        self._n_y = max(1, n_y)

    def draw(self, renderer) -> None:
        if self._min_cell_pt > 0 and self.axes is not None and self.figure is not None:
            bbox = self.axes.get_window_extent(renderer)
            cell_width = bbox.width * 72 / self.figure.dpi / self._n_x
            cell_height = bbox.height * 72 / self.figure.dpi / self._n_y
            if cell_width < self._min_cell_pt or cell_height < self._min_cell_pt:
                return
        super().draw(renderer)


class geom_contrast_text(geom_text):
    """Draw vector cell text only when legible, with fill-aware contrast."""

    DEFAULT_AES = {**geom_text.DEFAULT_AES, "fill": "white"}
    DEFAULT_PARAMS = {
        **geom_text.DEFAULT_PARAMS,
        "min_cell_pt": 12.0,
        "format_string": "{}",
    }

    def draw_panel(self, data, panel_params, coord, ax) -> None:
        # Text inherited from a continuous fill would otherwise be split into
        # one group per colour, making every group look like a one-cell panel.
        self.draw_group(data, panel_params, coord, ax, self.params)

    @staticmethod
    def draw_group(
        data: pd.DataFrame,
        panel_params,
        coord,
        ax,
        params: dict[str, Any],
    ) -> None:
        data = coord.transform(data, panel_params)
        n_x = int(data["x"].nunique(dropna=True))
        n_y = int(data["y"].nunique(dropna=True))
        for _, row in data.iterrows():
            fill = row.get("fill", "white")
            colour = "white" if _relative_luminance(fill) < 0.179 else "black"
            raw_label = row["label"]
            try:
                label = params["format_string"].format(raw_label)
            except (TypeError, ValueError):
                label = params["format_string"].format(float(raw_label))
            text = _CellText(
                x=row["x"],
                y=row["y"],
                text=label,
                min_cell_pt=float(params["min_cell_pt"]),
                n_x=n_x,
                n_y=n_y,
                color=colour,
                alpha=row.get("alpha", 1),
                size=row.get("size", 6),
                family=row.get("family", None),
                fontstyle=row.get("fontstyle", "normal"),
                fontweight=row.get("fontweight", "normal"),
                rotation=row.get("angle", 0),
                linespacing=row.get("lineheight", 1.2),
                ha=row.get("ha", "center"),
                va=row.get("va", "center"),
                zorder=params["zorder"],
                clip_on=True,
                transform=ax.transData,
            )
            ax._add_text(text)


def annotation_threshold(value: bool | str) -> float | None:
    """Translate public annotation modes into a rendered cell-size threshold."""
    if value is False:
        return None
    if value is True or value == "force":
        return 0.0
    if value == "auto":
        return 12.0
    raise ValueError("annotate must be False, True, 'auto', or 'force'.")


__all__ = ["geom_contrast_text", "annotation_threshold"]
