import astropy.units as u
import matplotlib.pyplot as plt
import sunpy.visualization.colormaps as cm  # NOQA
from mpl_animators import ArrayAnimatorWCS
from ndcube.visualization.mpl_sequence_plotter import MatplotlibSequencePlotter, SequenceAnimator

__all__ = ["IRISSequencePlotter", "IRISSequenceAnimator"]


def _set_axis_colors(ax):
    if isinstance(ax, ArrayAnimatorWCS):
        ax = ax.axes
    if len(ax.coords._as_table()) == 1:
        wave = ax.coords[0]
        wave.set_format_unit(u.nm)
        wave.set_major_formatter("x.x")
        return
    if len(ax.coords._as_table()) == 2:
        lon, lat = ax.coords
    elif len(ax.coords._as_table()) == 3:
        wave, lat, lon = ax.coords
        wave.set_format_unit(u.nm)
        wave.set_major_formatter("x.x")
    else:
        raise ValueError(f"Too many axes: {len(ax.coords._as_table())}")
    lon.set_ticklabel_position("all")
    lat.set_ticklabel_position("all")
    lon.set_axislabel(ax.get_xlabel(), color="black")
    lon.set_ticklabel("black")
    lat.set_axislabel(ax.get_ylabel(), color="red")
    lat.set_ticklabel("red")


class IRISSequencePlotter(MatplotlibSequencePlotter):
    def plot(self, sequence_axis_coords=None, sequence_axis_unit=None, **kwargs):
        cmap = kwargs.get("cmap")
        if not cmap:
            try:
                cmap = plt.get_cmap(name="irissji{}".format(int(self.meta["TWAVE1"])))
            except Exception:
                cmap = "viridis"
        kwargs["cmap"] = cmap
        return IRISSequenceAnimator(self, sequence_axis_coords, sequence_axis_unit, **kwargs)

    def animate(self, sequence_axis_coords=None, sequence_axis_unit=None, **kwargs):
        cmap = kwargs.get("cmap")
        if not cmap:
            try:
                cmap = plt.get_cmap(name="irissji{}".format(int(self.meta["TWAVE1"])))
            except Exception:
                cmap = "viridis"
        kwargs["cmap"] = cmap
        return IRISSequenceAnimator(self, sequence_axis_coords, sequence_axis_unit, **kwargs)


class IRISSequenceAnimator(SequenceAnimator):
    def plot_start_image_2d(self, ax):
        im = super().plot_start_image_2d(ax)
        _set_axis_colors(ax)
        return im

    def update_plot_2d(self, val, im, slider):
        super().update_plot_2d(val, im, slider)
        _set_axis_colors(self.axes)
