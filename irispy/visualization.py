import matplotlib.pyplot as plt
import sunpy.visualization.colormaps as cm  # NOQA
from ndcube.visualization.mpl_sequence_plotter import MatplotlibSequencePlotter, SequenceAnimator

__all__ = ["IRISSequencePlotter", "IRISSequenceAnimator"]


def _set_axis_colors(ax):
    lon, lat = ax.coords
    lon.set_ticklabel_position("all")
    lat.set_ticklabel_position("all")
    lon.set_axislabel(ax.get_xlabel(), color="tab:blue")
    lon.set_ticklabel("tab:blue")
    lat.set_axislabel(ax.get_ylabel(), color="tab:green")
    lat.set_ticklabel("tab:green")


class IRISSequencePlotter(MatplotlibSequencePlotter):
    def plot(self, sequence_axis_coords=None, sequence_axis_unit=None, **kwargs):
        cmap = kwargs.get("cmap")
        if not cmap:
            cmap = plt.get_cmap(name="irissji{}".format(int(self.meta["TWAVE1"])))
        kwargs["cmap"] = cmap
        return IRISSequenceAnimator(self, sequence_axis_coords, sequence_axis_unit, **kwargs)

    def animate(self, sequence_axis_coords=None, sequence_axis_unit=None, **kwargs):
        cmap = kwargs.get("cmap")
        if not cmap:
            cmap = plt.get_cmap(name="irissji{}".format(int(self.meta["TWAVE1"])))
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