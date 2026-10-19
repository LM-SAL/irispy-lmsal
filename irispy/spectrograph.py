import textwrap
from copy import deepcopy
from numbers import Integral

import numpy as np

import astropy.units as u
import gwcs
from astropy.time import Time
from astropy.wcs import WCS
from astropy.wcs.wcsapi import SlicedLowLevelWCS

from ndcube import NDCollection
from ndcube.utils.wcs import pixel_axis_to_physical_types
from ndcube.visualization import PlotterDescriptor
from ndcube.wcs.tools import unwrap_wcs_to_fitswcs
from ndcube.wcs.wrappers import ResampledLowLevelWCS
from sunpy import log as logger
from sunraster import SpectrogramCube as SpecCube
from sunraster.spectrogram import _calculate_exposure_time_correction, _uncalculate_exposure_time_correction

from irispy._spectrograph_wcs import _raster_crop_bounds
from irispy._wcs import _celestial_frame_from_cube, _ResolveNegativeIndicesMixin
from irispy.utils.constants import SLIT_WIDTH
from irispy.utils.cosmic_rays import remove_cosmic_rays
from irispy.visualization import SpectrogramPlotter

__all__ = ["MosaicCube", "RasterCollection", "SpectrogramCube"]


class SpectrogramCube(_ResolveNegativeIndicesMixin, SpecCube):
    """
    Class representing spectrogram data described by a single WCS.

    A raster window is exposed as one cube, whether it comes from a single file,
    a combined multi-file raster, or a sit-and-stare observation.

    Parameters
    ----------
    data: `numpy.ndarray`
        The array holding the actual data in this object.
    wcs: `astropy.wcs.WCS` or ``gwcs.WCS``
        The WCS object containing the axes' information
    unit : `astropy.units.Unit` or `str`, optional
        Unit for the dataset. Strings that can be converted to a Unit are allowed.
        Defaults to None, which takes the unit of ``data`` if it is a `~astropy.units.Quantity`.
    meta : `dict` object, optional
        Additional meta information about the dataset, such as the
        `~irispy.meta.SGMeta` of the observation the data comes from.
        Defaults to None.
    uncertainty : any type, optional
        Uncertainty in the dataset. Should have an attribute uncertainty_type
        that defines what kind of uncertainty is stored, for example "std"
        for standard deviation or "var" for variance. A metaclass defining
        such an interface is NDUncertainty - but isn't mandatory. If the uncertainty
        has no such attribute the uncertainty is stored as UnknownUncertainty.
        Defaults to None.
    mask : any type, optional
        Mask for the dataset. Masks should follow the numpy convention
        that valid data points are marked by False and invalid ones with True.
        Defaults to None.
    copy : `bool`, optional
        Indicates whether to save the arguments as copy. True copies every attribute
        before saving it while False tries to save every parameter as reference.
        Note however that it is not always possible to save the input as reference.
        Default is False.

    Examples
    --------
    The array axes are raster step or time, slit position and wavelength; the exposure
    times are an extra coordinate:

    >>> from irispy.data.test import get_test_filepath
    >>> from irispy.io import read_files
    >>> filename = get_test_filepath(
    ...     "sns/iris_l2_20210905_001833_3620258102_raster_t000_r00000_test.fits"
    ... )
    >>> cube = read_files(filename, spectral_windows="Si IV 1403")["Si IV 1403"]
    >>> cube.shape, cube.unit, cube.wavelength_axis
    ((187, 40, 29), Unit("DN_IRIS_FUV"), 2)
    >>> cube.meta.spectral_window, cube.meta.rest_wavelength
    ('Si IV 1403', <Quantity 140.277 nm>)
    >>> cube.spectral_dispersion.to("pm")
    <Quantity 25.44 pm>
    >>> cube[:, 10:20, 5:10].shape
    (187, 10, 5)
    """

    plotter = PlotterDescriptor(default_type=SpectrogramPlotter)

    def __init__(self, data, wcs, uncertainty=None, unit=None, meta=None, *, mask=None, copy=False, **kwargs) -> None:
        super().__init__(data, wcs, unit=unit, uncertainty=uncertainty, mask=mask, meta=meta, copy=copy, **kwargs)

    def __repr__(self) -> str:
        return f"{object.__repr__(self)}\n{self!s}"

    def _time_bounds(self):
        if "time" in self.global_coords:
            return self.global_coords["time"].min().isot, self.global_coords["time"].max().isot
        try:
            extra_coord_time = self.axis_world_coords("time", wcs=self.extra_coords)
            if extra_coord_time:
                return extra_coord_time[0].min().isot, extra_coord_time[0].max().isot
        except ValueError as e:
            logger.debug(f"Unable to determine time bounds for SpectrogramCube string representation: {e}")
        try:
            return self.time.min().isot, self.time.max().isot
        except (ValueError, AttributeError) as e:
            logger.debug(f"Unable to determine time bounds for SpectrogramCube string representation: {e}")
            return "Unknown", "Unknown"

    def __str__(self) -> str:
        instance_start, instance_end = self._time_bounds()
        return textwrap.dedent(
            f"""
            SpectrogramCube
            ---------------
            Obs ID:             {self.meta.get("OBSID")}
            Obs Description:    {self.meta.get("OBS_DESC")}
            Obs Date:           {instance_start} -- {instance_end}
            Data shape:         {self.shape}
            Axis Types:         {self._wcs_array_axis_physical_types}
            Roll:               {self.meta.get("SAT_ROT")}
            """,
        )

    remove_cosmic_rays = remove_cosmic_rays

    celestial_frame = property(_celestial_frame_from_cube)
    _get_crop_bounds = _raster_crop_bounds

    def axis_world_coords(self, *axes, **kwargs):
        # TODO: remove once ndcube returns C-ordered world coordinates or astropy's Time string
        # formats respect memory order. ndcube transposes N-D coordinates (Fortran order), and
        # astropy's isot/iso/datetime64 then emit them out of order, scrambling a 4D cube's times.
        # astropy main emits them in order since astropy/astropy#19942, which is not released yet.
        return tuple(
            coord.copy() if isinstance(coord, Time) else coord for coord in super().axis_world_coords(*axes, **kwargs)
        )

    @property
    def fits_wcs(self):
        """
        The plain FITS WCS built from the window header, the WCS of a cube built on a
        FITS WCS (such as a mosaic), or `None` when no single FITS WCS describes this
        cube (for example a combined multi-file cube).
        """
        root, low = self.meta.get("fits_wcs") if hasattr(self.meta, "get") else None, self.wcs.low_level_wcs
        item = [slice(None)] * low.pixel_n_dim
        if isinstance(low, SlicedLowLevelWCS):
            low, item = low._wcs, list(low._slices_array)
        if isinstance(low, WCS):
            root = low
        elif root is None or not isinstance(low, gwcs.WCS):
            return None
        if isinstance(root, list):
            if not isinstance(item[0], Integral):
                return None
            root, item = root[item[0]], item[1:]
        return root if all(i == slice(None) for i in item) else root.slice(tuple(item), numpy_order=True)

    @property
    def _separate_raster_axis(self):
        return "raster_scan" in self.wcs.low_level_wcs.world_axis_names

    def raster_slice(self, index):
        """
        Return the subcube corresponding to one original raster.
        """
        if not isinstance(index, Integral):
            msg = "Raster index must be an integer."
            raise TypeError(msg)

        n_rasters = self.shape[0] if self._separate_raster_axis else 1
        if not -n_rasters <= index < n_rasters:
            msg = "Raster index out of range."
            raise IndexError(msg)

        return self[index] if self._separate_raster_axis else self

    def split_rasters(self):
        """
        Split the cube into per-raster subcubes.
        """
        return tuple(self[i] for i in range(self.shape[0])) if self._separate_raster_axis else (self,)

    def apply_exposure_time_correction(self, undo=False, force=False):  # NOQA: FBT002 (sunraster signature)
        exposure_time = self.exposure_time.to_value(u.s)
        if np.ndim(exposure_time) < 2:
            return super().apply_exposure_time_correction(undo=undo, force=force)
        # sunraster expects one exposure axis; a combined cube has two (scan and step).
        axes = self._get_axis_coord_index(self._exposure_time_name, self._exposure_time_loc)
        item = tuple(slice(None) if axis in axes else np.newaxis for axis in range(self.data.ndim))
        correct = _uncalculate_exposure_time_correction if undo else _calculate_exposure_time_correction
        new_cube = deepcopy(self)
        new_cube._data, new_cube._uncertainty, new_cube._unit = correct(
            self.data, self.uncertainty, self.unit, exposure_time[item], force=force
        )
        return new_cube

    @property
    def _wcs_array_axis_physical_types(self):
        """
        `~ndcube.NDCube.array_axis_physical_types` from the main WCS alone.

        The combined WCS, with the extra coordinates, raises on a length-1 step axis
        (ndcube: mismatched pixel bounds), such as ``cube[0:1]`` or a one-exposure file.
        """
        wcs = self.wcs.low_level_wcs
        return [tuple(pixel_axis_to_physical_types(axis, wcs).tolist()) for axis in reversed(range(wcs.pixel_n_dim))]

    @property
    def _fits_wcsprm(self):
        """
        Raw FITS WCS keywords (``Wcsprm``) of the native pixels of this cube.

        Slicing and rebinning keep the native pixel scales, so these come from the
        unsliced WCS: ``self.wcs`` for a cube built on a FITS WCS, otherwise the
        per-file ``fits_wcs`` in ``meta``.
        """
        root = self.wcs.low_level_wcs
        while isinstance(root, (SlicedLowLevelWCS, ResampledLowLevelWCS)):
            root = root._wcs
        fits_wcs = root if hasattr(root, "wcs") else self.meta.get("fits_wcs") if hasattr(self.meta, "get") else None
        fits_wcs = fits_wcs[0] if isinstance(fits_wcs, list) else fits_wcs
        if fits_wcs is None:
            msg = "This cube has no FITS WCS to take its pixel scales from."
            raise ValueError(msg)
        return fits_wcs.wcs

    @property
    def spectral_dispersion(self):
        """
        Spectral dispersion per pixel along the wavelength axis.
        """
        wcs = self._fits_wcsprm
        mask = np.array([ctype == "WAVE" for ctype in wcs.ctype])
        if not mask.any():
            msg = "Cannot determine spectral axis (no WAVE ctype in WCS) for spectral_dispersion"
            raise ValueError(msg)
        idx = np.argmax(mask)
        return wcs.cdelt[idx] * u.Unit(wcs.cunit[idx])

    @property
    def solid_angle(self):
        """
        Solid angle per spatial pixel (slit width x spatial pixel scale).
        """
        wcs = self._fits_wcsprm
        mask = np.array(["HPLT" in ctype for ctype in wcs.ctype])
        if not mask.any():
            msg = "Cannot determine latitude axis (no HPLT ctype in WCS) for solid_angle computation"
            raise ValueError(msg)
        lat_idx = np.argmax(mask)
        return wcs.cdelt[lat_idx] * u.Unit(wcs.cunit[lat_idx]) * SLIT_WIDTH

    @property
    def wavelength_axis(self):
        """
        Index of the spectral (wavelength) axis.
        """
        try:
            return next(
                axis
                for axis, physical_types in enumerate(self._wcs_array_axis_physical_types)
                if physical_types and "em.wl" in physical_types
            )
        except StopIteration:
            msg = "Could not identify a spectral wavelength axis on the cube"
            raise ValueError(msg) from None


def _wavelength_indices(wavelengths, wavelength_range, *, allow_empty=False):
    """
    Indices of the ``wavelengths`` between two wavelengths, inclusive.

    The bounds are widened by 1e-6 Å, so that a wavelength typed from the grid is kept
    despite the rounding of the WCS unit conversion. With ``allow_empty=True``, a range
    without samples returns an empty array.
    """
    wavelength_range = u.Quantity(wavelength_range, u.AA)
    if wavelength_range.shape != (2,):
        msg = f"Expected two wavelengths, got {wavelength_range}"
        raise ValueError(msg)
    low, high = wavelength_range + [-1e-6, 1e-6] * u.AA
    indices = np.flatnonzero((wavelengths >= low) & (wavelengths <= high))
    if not indices.size and not allow_empty:
        msg = f"No wavelengths between {wavelength_range[0]} and {wavelength_range[1]}"
        raise ValueError(msg)
    return indices


class MosaicCube(SpectrogramCube):
    """
    An IRIS full-disk mosaic of one spectral window, read by `irispy.io.read_mosaic`.

    The array axes are (wavelength, solar Y, solar X), the reverse of a Level 2 raster,
    so the raster analysis helpers reject it; use `to_map` and
    `~irispy.utils.moments.calculate_moments`.
    """

    def __str__(self) -> str:
        return textwrap.dedent(
            f"""
            MosaicCube
            ----------
            Obs ID:          {self.meta.get("OBSID")}
            Rest wavelength: {self.meta.get("LAMREF")} Angstrom
            Obs Date:        {self.meta.get("DATE_OBS")} -- {self.meta.get("DATE_END")}
            Data shape:      {self.shape}
            Axis Types:      {self.array_axis_physical_types}
            """,
        )

    def to_map(self, wavelength):
        """
        Return a `sunpy.map.Map` at one wavelength or averaged over a range.

        The color map is that of the nearest slit-jaw passband. Masked values are left
        out of the average.

        Parameters
        ----------
        wavelength : `astropy.units.Quantity`
            One wavelength, for the nearest mosaic wavelength, or two, for the mean
            between them (inclusive).

        Returns
        -------
        `sunpy.map.GenericMap`
        """
        from sunpy.map import Map  # NOQA: PLC0415

        wavelength = u.Quantity(wavelength, u.AA)
        axis = self.wavelength_axis
        wavelengths = u.Quantity(self.axis_world_coords(axis)[0], u.AA)
        if wavelength.isscalar:
            indices = [np.argmin(np.abs(wavelengths - wavelength))]
        else:
            indices = _wavelength_indices(wavelengths, wavelength)
        data = np.take(self.data, indices, axis=axis)
        if self.mask is not None:
            data = np.where(np.take(self.mask, indices, axis=axis), np.nan, data)
        # A masked mean leaves positions with no unmasked value NaN, without a warning.
        data = np.ma.masked_invalid(data).mean(axis=axis).filled(np.nan)
        header = unwrap_wcs_to_fitswcs(self.wcs)[0].celestial.to_header()
        mean_wavelength = wavelengths[indices].mean()
        header.update(
            {
                "TELESCOP": "IRIS",
                "OBSRVTRY": "IRIS",
                "INSTRUME": "SPEC",  # the Level 2 spectrograph value
                "EXPTIME": self.meta.get("EXPTIME"),
                "DATE-END": self.meta.get("DATE_END"),
                "WAVELNTH": mean_wavelength.to_value(u.AA),
                "WAVEUNIT": "Angstrom",
                "BUNIT": "DN",
            }
        )
        sunpy_map = Map(data, header)
        passband = min((1330, 1400, 2796, 2832), key=lambda band: abs(band - mean_wavelength.to_value(u.AA)))
        sunpy_map.plot_settings["cmap"] = f"irissji{passband}"
        return sunpy_map


class RasterCollection(NDCollection):
    """
    Subclass of NDCollection for raster spectral windows keyed by window name.

    Each value is a `SpectrogramCube`.
    """

    def __str__(self) -> str:
        return textwrap.dedent(
            f"""
            Raster Collection
            -----------------
            Spectral Windows (cube keys): {tuple(self.keys())}
            Number of Cubes: {len(self)}
            Aligned dimensions: {self.aligned_dimensions}
            Aligned physical types: {self.aligned_axis_physical_types}
            """,
        )

    @property
    def aligned_axis_physical_types(self):
        # NDCollection builds each tuple from a set, whose order changes from run to run.
        # TODO: remove once irispy requires an ndcube with sunpy/ndcube#983.
        types = super().aligned_axis_physical_types
        return None if types is None else [tuple(sorted(axis_types)) for axis_types in types]
