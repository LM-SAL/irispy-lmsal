import astropy.units as u
import sunpy.visualization.colormaps as cm  # NOQA: F401
from mpl_animators import ArrayAnimatorWCS
from ndcube.visualization.mpl_plotter import MatplotlibPlotter

__all__ = ["_set_axis_colors", "Plotter", "CustomArrayAnimatorWCS"]


def _set_axis_colors(ax):
    """
    Private method to set the axis colors and labels for IRIS SJI and Raster
    data.
    """
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
        _set_axis_colors(self.axes)


class Plotter(MatplotlibPlotter):
    def plot(self, *args, **kwargs):
        return super().plot(*args, **kwargs)

    def _animate_cube(self, wcs, plot_axes=None, axes_coordinates=None, axes_units=None, data_unit=None, **kwargs):
        # Derive inputs for animation object and instantiate.
        data, wcs, plot_axes, coord_params = self._prep_animate_args(wcs, plot_axes, axes_units, data_unit)
        ax = CustomArrayAnimatorWCS(data, wcs, plot_axes, coord_params=coord_params, **kwargs)
        # We need to modify the visible axes after the axes object has been created.
        # This call affects only the initial draw
        self._apply_axes_coordinates(ax.axes, axes_coordinates)
        # This changes the parameters for future iterations
        for hidden in self._not_visible_coords(ax.axes, axes_coordinates):
            if hidden in ax.coord_params:
                param = ax.coord_params[hidden]
            else:
                param = {}
            param["ticks"] = False
            ax.coord_params[hidden] = param
        return ax
