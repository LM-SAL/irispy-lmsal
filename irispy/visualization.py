import matplotlib as mpl
import numpy as np
from mpl_animators import ArrayAnimatorWCS

import astropy.units as u
from astropy.visualization import AsymmetricPercentileInterval
from astropy.wcs.wcsapi import HighLevelWCSWrapper

import sunpy.visualization.colormaps as cm  # NOQA: F401
from ndcube.visualization.mpl_plotter import MatplotlibPlotter
from ndcube.visualization.mpl_sequence_plotter import MatplotlibSequencePlotter, SequenceAnimator

__all__ = []


LAT_LABELS = [
    "custom:pos.helioprojective.lat",
    "helioprojective latitude",
    "hplt-tan",
    "hplt",
    "lat",
    "latitude",
]
LON_LABELS = [
    "custom:pos.helioprojective.lon",
    "helioprojective longitude",
    "hpln-tan",
    "hpln",
    "lon",
    "longitude",
]
TIME_LABELS = ["time (utc)", "time"]
WAVELENGTH_LABELS = ["wavelength", "wave", "em.wl"]
SCAN_LABELS = set(LON_LABELS + LAT_LABELS)
LON_AXIS_LABEL = "Helioprojective Longitude [arcsec]"
LAT_AXIS_LABEL = "Helioprojective Latitude [arcsec]"


def _shorten_slider_label(label):
    label_parts = {label_part.strip().lower() for label_part in str(label).split(" / ")}
    if label_parts.intersection(WAVELENGTH_LABELS):
        return "Wavelength"
    if label_parts.intersection(TIME_LABELS):
        return "Time"
    if label_parts and label_parts.issubset(SCAN_LABELS):
        return "Scan"
    return label


def set_axis_properties(ax):
    """
    Set the axis colors and labels for IRIS SJI and Raster data.
    """
    if isinstance(ax, ArrayAnimatorWCS):
        ax = ax.axes
    for axis in ax.coords:
        if axis.default_label.lower() in WAVELENGTH_LABELS:
            axis.set_format_unit(u.nm)
            axis.set_major_formatter("x.xx")
            # A narrow window drawn with equal aspect has no room to label every tick
            axis.set_ticklabel(exclude_overlapping=True)
            axis.set_axislabel("Wavelength [$\\mathrm{nm}$]")
        elif axis.default_label.lower() in LAT_LABELS:
            _set_axis_properties(axis, LAT_AXIS_LABEL, "red")
        elif axis.default_label.lower() in LON_LABELS:
            _set_axis_properties(axis, LON_AXIS_LABEL, "black")


def _set_axis_properties(axis, label, color):
    """
    Set the axis colors and labels for IRIS SJI and Raster data.

    Parameters
    ----------
    axis : `~astropy.visualization.wcsaxes.core.WCSAxes`
        The axis to set the colors and labels for.
    label : str
        The label to use for the axis.
    color : str
        The color to use for the axis label.
    """
    axis.set_ticklabel(color, fontsize=8)
    axis.set_axislabel(label, color=color, fontsize=8)


def _map_aspect(cube, wcs, plot_axes):
    """
    The ``aspect`` that draws a 2D map with its true proportions, `None` for other
    images.

    ``plot_axes`` is in pixel order, as `ndcube` passes it to ``_plot_2D_cube``.
    """
    spatial = set(wcs.world_axis_physical_types) == {"custom:pos.helioprojective.lon", "custom:pos.helioprojective.lat"}
    # The step axis of a sit-and-stare is time, along which the slit only follows solar rotation; without
    # metadata a cube cannot tell one from a raster.
    if not spatial or not cube.meta or cube.meta.get("NRASTERP") == 1:
        return None
    nx, ny = cube.data.shape[::-1]
    wcs = HighLevelWCSWrapper(wcs)
    scale = []
    # The mean pixel spacing across the middle of the map, which pointing jitter barely changes. Steps
    # without an exposure have no pointing, so measure between the finite ends of each line.
    for line in (
        wcs.pixel_to_world(np.arange(nx), np.full(nx, ny // 2)),
        wcs.pixel_to_world(np.full(ny, nx // 2), np.arange(ny)),
    ):
        finite = np.flatnonzero(np.isfinite(line.Tx))
        if finite.size < 2:
            return None
        scale.append(line[finite[0]].separation(line[finite[-1]]) / (finite[-1] - finite[0]))
    return (scale[plot_axes.index("y")] / scale[plot_axes.index("x")]).to_value(u.one)


class Plot2DMixin:
    def update_plot_2d(self, val, im, slider):
        super().update_plot_2d(val, im, slider)
        set_axis_properties(self.axes)


class IRISArrayAnimatorWCS(Plot2DMixin, ArrayAnimatorWCS):
    def _compute_slider_labels_from_wcs(self, slices):
        return [_shorten_slider_label(label) for label in super()._compute_slider_labels_from_wcs(slices)]


class IRISSequenceAnimator(Plot2DMixin, SequenceAnimator):
    def __init__(self, *args, **kwargs):
        kwargs.setdefault("slider_labels", ["Raster step", "Scan number"])
        super().__init__(*args, **kwargs)


def _registered_colormap(name):
    """
    ``name`` when matplotlib knows that colormap, viridis otherwise.
    """
    return name if name in mpl.colormaps else "viridis"


def _make_map(data, header, unit, mask):
    """
    A `sunpy.map.Map` of ``data`` and ``mask`` with a copy of ``header`` that gives
    ``unit`` as its BUNIT.
    """
    from sunpy.map import Map  # NOQA: PLC0415

    # The IRIS DN units, such as DN_IRIS_SJI_1330, are not FITS units; sunpy reads DN
    bases = [u.DN if base.name.startswith("DN_IRIS") else base for base in unit.bases]
    header = header.copy()
    header["bunit"] = u.CompositeUnit(unit.scale, bases, unit.powers).to_string()
    return Map(data, header, mask=None if mask is None else np.asarray(mask))


def _sji_cmap_name(meta):
    """
    The sunpy color table of a slit-jaw or AIA cube's passband, viridis without one.
    """
    # Hinode/SOT cubes have TWAVE1 = 0, and a cube built by hand may have none
    twave = meta.get("TWAVE1")
    if not twave:
        return "viridis"
    # sunpy registers the AIA tables as sdoaia<channel>. Its AIA maps build that name from TELESCOP,
    # which is blank in the AIA files supplied with IRIS, and ask for the unregistered aia<channel>.
    prefix = "sdoaia" if meta.get("INSTRUME", "").startswith("AIA") else "irissji"
    return _registered_colormap(f"{prefix}{int(twave)}")


class _IRISPlotMixin:
    """
    Shared IRIS plot behavior for cube and sequence plotters.
    """

    def _default_cmap_name(self):
        # sunpy registers the IRIS spectrograph color tables as irissjiFUV and irissjiNUV
        return _registered_colormap(f"irissji{getattr(self._ndcube.meta, 'detector_band', '')}")

    def plot(self, *args, **kwargs):
        """
        Plot the cube with IRIS defaults.

        If ``cmap`` is falsey, the color table of the data is used: the IRIS passband
        table for slit-jaw images, the sunpy ``sdoaia`` table of the channel for AIA
        cubes and the IRIS FUV or NUV table for spectrograms, with viridis for
        Hinode/SOT cubes. For one-dimensional cubes ``cmap`` is dropped entirely, even
        when given. IRIS axis styling is applied to the result; all other arguments are
        passed to the parent plotter's ``plot``.

        A 2D map, such as a cut of a scanning raster at one wavelength, a moment map or
        a Hinode/SOT-SP scan, is drawn with its true proportions unless ``aspect`` is
        given. Maps of a sit-and-stare keep square pixels.

        By default the colors run linearly from the minimum to the maximum of the image,
        or of the first frame of an animation, which the other frames keep, so a few
        bright pixels, such as cosmic-ray hits, can leave little contrast.

        ``clip_interval=(1, 99.5) * u.percent`` sets them to those percentiles of the
        image, or of each frame.
        """
        kwargs["cmap"] = kwargs.get("cmap") or self._default_cmap_name()
        if len(self._ndcube.shape) == 1:
            kwargs.pop("cmap")
        ax = super().plot(*args, **kwargs)
        set_axis_properties(ax)
        return ax


class IRISPlotter(_IRISPlotMixin, MatplotlibPlotter):
    def _plot_2D_cube(  # NOQA: N802
        self, wcs, axes=None, plot_axes=None, axes_coordinates=None, axes_units=None, data_unit=None, **kwargs
    ):
        if "aspect" not in kwargs:
            kwargs["aspect"] = _map_aspect(self._ndcube, wcs, plot_axes)
        # ndcube passes clip_interval on to imshow, which rejects it; remove this once the minimum
        # ndcube applies it to 2D cubes, as its animations do.
        if (clip_interval := kwargs.pop("clip_interval", None)) is not None:
            if len(clip_interval) != 2:
                msg = "A range of 2 values must be specified for clip_interval."
                raise ValueError(msg)
            data = self._ndcube.data if data_unit is None else self._ndcube.quantity.to_value(data_unit)
            interval = AsymmetricPercentileInterval(*clip_interval.to_value(u.percent))
            kwargs["vmin"], kwargs["vmax"] = interval.get_limits(np.ma.masked_array(data, self._ndcube.mask))
        return super()._plot_2D_cube(wcs, axes, plot_axes, axes_coordinates, axes_units, data_unit, **kwargs)

    def _animate_cube(
        self,
        wcs,
        plot_axes=None,
        axes_coordinates=None,
        axes_units=None,
        data_unit=None,
        **kwargs,
    ):
        # This is a copy of the super method, but with the replacement
        # of the ArrayAnimatorWCS with IRISArrayAnimatorWCS.
        data, wcs, plot_axes, coord_params = self._prep_animate_args(wcs, plot_axes, axes_units, data_unit)
        ax = IRISArrayAnimatorWCS(data, wcs, plot_axes, coord_params=coord_params, **kwargs)
        # We need to modify the visible axes after the axes object has been created.
        # This call affects only the initial draw
        self._apply_axes_coordinates(ax.axes, axes_coordinates)
        # This changes the parameters for future iterations
        for hidden in self._not_visible_coords(ax.axes, axes_coordinates):
            param = ax.coord_params.get(hidden, {})
            param["ticks"] = False
            ax.coord_params[hidden] = param
        return ax


class SJIPlotter(IRISPlotter):
    def _default_cmap_name(self):
        return _sji_cmap_name(self._ndcube.meta)


class IRISSequencePlotter(_IRISPlotMixin, MatplotlibSequencePlotter):
    def animate(self, sequence_axis_coords=None, sequence_axis_unit=None, **kwargs):
        """
        Animate the `~ndcube.NDCubeSequence` with the sequence axis as a slider.

        Keyword arguments are passed to
        `ndcube.visualization.mpl_plotter.MatplotlibPlotter.plot` and therefore only
        apply to cube axes, not the sequence axis.
        See that method's docstring for definition of keyword arguments.

        Parameters
        ----------
        sequence_axis_coords: `str` optional
            The name of the coordinate in `~ndcube.NDCubeSequence.sequence_axis_coords`
            to be used as the slider pixel values.
            If None, array indices will be used.

        sequence_axis_unit: `astropy.units.Unit` or `str`, optional
            The unit in which the sequence_axis_coordinates should be displayed.
            If None, the default unit will be used.
        """
        return IRISSequenceAnimator(
            self._ndcube, sequence_axis_coords=sequence_axis_coords, sequence_axis_unit=sequence_axis_unit, **kwargs
        )
