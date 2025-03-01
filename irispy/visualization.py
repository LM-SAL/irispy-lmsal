from mpl_animators import ArrayAnimatorWCS

import astropy.units as u

import sunpy.visualization.colormaps as cm  # NOQA: F401
from ndcube.visualization.mpl_plotter import MatplotlibPlotter

__all__ = ["CustomArrayAnimatorWCS", "Plotter", "set_axis_properties"]


LAT_LABELS = [
    "custom:pos.helioprojective.lat",
    "hplt-tan",
    "hplt",
    "lat",
    "latitude",
]
LON_LABELS = [
    "custom:pos.helioprojective.lon",
    "hpln-tan",
    "hpln",
    "lon",
    "longitude",
]
WAVELENGTH_LABELS = ["wavelength", "wave", "em.wl"]


def set_axis_properties(ax):
    """
    Set the axis colors and labels for IRIS SJI and Raster data.
    """
    if isinstance(ax, ArrayAnimatorWCS):
        ax = ax.axes
    for axis in ax.coords:
        if axis.default_label.lower() in WAVELENGTH_LABELS:
            axis.set_format_unit(u.nm)
            axis.set_major_formatter("x.x")
            axis.set_axislabel("Wavelength [$\\mathrm{nm}$]")
        elif axis.default_label.lower() in LAT_LABELS:
            _set_axis_properties(axis, "red")
        elif axis.default_label.lower() in LON_LABELS:
            _set_axis_properties(axis, "black")


def _set_axis_properties(axis, color):
    """
    Set the axis colors and labels for IRIS SJI and Raster data.

    Parameters
    ----------
    axis : `~astropy.visualization.wcsaxes.core.WCSAxes`
        The axis to set the colors and labels for.
    color : str
        The color to use for the axis label.
    """
    axis.set_ticklabel(color, fontsize=8)
    axis.set_axislabel(axis.default_label, color=color, fontsize=8)


class CustomArrayAnimatorWCS(ArrayAnimatorWCS):
    def update_plot_2d(self, val, im, slider):
        """
        Update the image plot.
        """
        self.axes.reset_wcs(wcs=self.wcs, slices=self.slices_wcsaxes)
        im.set_array(self.data_transposed)
        if self.clip_interval is not None:
            vmin, vmax = self._get_2d_plot_limits()
            im.set_clim(vmin, vmax)
        slider.cval = val
        set_axis_properties(self.axes)


class Plotter(MatplotlibPlotter):
    def plot(self, *args, **kwargs):
        return super().plot(*args, **kwargs)

    def _animate_cube(
        self,
        wcs,
        plot_axes=None,
        axes_coordinates=None,  # NOQA: ARG002 - Unused but passed in via ModestImage.set
        axes_units=None,
        data_unit=None,
        **kwargs,
    ):
        data, wcs, plot_axes, coord_params = self._prep_animate_args(wcs, plot_axes, axes_units, data_unit)
        return CustomArrayAnimatorWCS(data, wcs, plot_axes, coord_params=coord_params, **kwargs)
