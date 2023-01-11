import astropy.units as u
import matplotlib.pyplot as plt
import sunpy.visualization.colormaps as cm  # NOQA
from mpl_animators import ArrayAnimatorWCS
from ndcube.visualization.mpl_sequence_plotter import MatplotlibSequencePlotter
from ndcube.visualization.mpl_sequence_plotter import SequenceAnimator as SeqAnim

__all__ = ["SequencePlotter", "SequenceAnimator"]


def _set_axis_colors(ax):
    if isinstance(ax, ArrayAnimatorWCS):
        ax = ax.axes
    for axis in ax.coords:
        if axis.default_label.lower() in ["wavelength", "wave", "em.wl"]:
            axis.set_format_unit(u.nm)
            axis.set_major_formatter("x.x")
            axis.set_axislabel("em.wl [$\\mathrm{nm}$]")
        elif axis.default_label.lower() in ["latitude", "lat", "custom:pos.helioprojective.lat"]:
            axis.set_axislabel(axis.default_label, color="red")
            axis.set_ticklabel_position("all")
            axis.set_ticklabel("red")
        elif axis.default_label.lower() in ["longitude", "lon", "custom:pos.helioprojective.lon"]:
            axis.set_axislabel(axis.default_label, color="black")
            axis.set_ticklabel_position("all")
            axis.set_ticklabel("black")


class SequencePlotter(MatplotlibSequencePlotter):
    def plot(self, sequence_axis_coords=None, sequence_axis_unit=None, **kwargs):
        cmap = kwargs.get("cmap")
        if not cmap:
            try:
                cmap = plt.get_cmap(name="irissji{}".format(int(self.meta["TWAVE1"])))
            except Exception:
                cmap = "viridis"
        kwargs["cmap"] = cmap
        return SequenceAnimator(self, sequence_axis_coords, sequence_axis_unit, **kwargs)

    def animate(self, sequence_axis_coords=None, sequence_axis_unit=None, **kwargs):
        cmap = kwargs.get("cmap")
        if not cmap:
            try:
                cmap = plt.get_cmap(name="irissji{}".format(int(self.meta["TWAVE1"])))
            except Exception:
                cmap = "viridis"
        kwargs["cmap"] = cmap
        return SequenceAnimator(self, sequence_axis_coords, sequence_axis_unit, **kwargs)


class SequenceAnimator(SeqAnim):
    def plot_start_image_2d(self, ax):
        im = super().plot_start_image_2d(ax)
        _set_axis_colors(ax)
        return im

    def update_plot_2d(self, val, im, slider):
        super().update_plot_2d(val, im, slider)
        _set_axis_colors(self.axes)
