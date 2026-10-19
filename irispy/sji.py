import textwrap
from numbers import Integral
from functools import cached_property

import matplotlib as mpl
import numpy as np

import gwcs
from astropy.wcs import WCS
from astropy.wcs.wcsapi import SlicedLowLevelWCS

from ndcube.visualization import PlotterDescriptor
from sunpy.util import MetaDict
from sunraster import SpectrogramCube

from irispy._wcs import _celestial_frame_from_cube, _ResolveNegativeIndicesMixin
from irispy.utils import calculate_dust_mask
from irispy.utils._spectral import check_scaled
from irispy.utils.cosmic_rays import remove_cosmic_rays
from irispy.utils.dust import remove_dust as _remove_dust
from irispy.visualization import SJIPlotter

__all__ = ["AIACube", "SJICube", "SOTCube"]


def _slice_header(header, spatial):
    """
    The FITS header of one frame after the spatial slice ``spatial`` (numpy order) of
    its cube, in the file's own units.

    ``preserve_units`` keeps the header's arcsec keywords; without it ``to_header``
    rewrites the celestial keywords in degrees.
    """
    # ndcube spells a full axis slice(0, n, None)
    shape = (header["naxis2"], header["naxis1"])
    if all(isinstance(axis, slice) and axis.indices(n) == (0, n, 1) for axis, n in zip(spatial, shape, strict=True)):
        return header
    sliced = WCS(header, preserve_units=True).slice(spatial, numpy_order=True)
    header = header.copy()
    header.update(sliced.to_header())
    header["naxis1"], header["naxis2"] = sliced.pixel_shape
    return header


class SJICube(_ResolveNegativeIndicesMixin, SpectrogramCube):
    """
    Class representing SJI Image described by a single WCS.

    Parameters
    ----------
    data : `numpy.ndarray`
        The array holding the actual data in this object.
    wcs : `astropy.wcs.WCS`
        The WCS object containing the axes information
    unit : `astropy.units.Unit` or `str`
        Unit for the dataset.
        Strings that can be converted to a Unit are allowed.
    meta : `dict` object
        Additional meta information about the dataset.
    uncertainty : any type, optional
        Uncertainty in the dataset. Should have an attribute uncertainty_type
        that defines what kind of uncertainty is stored, for example "std"
        for standard deviation or "var" for variance. A metaclass defining
        such an interface is NDUncertainty - but isn't mandatory. If the
        uncertainty has no such attribute the uncertainty is stored as
        UnknownUncertainty.
        Defaults to None.
    mask : any type, optional
        Mask for the dataset. Masks should follow the numpy convention
        that valid data points are marked by False and invalid ones with True.
        Defaults to None.
    copy : `bool`, optional
        Indicates whether to save the arguments as copy. True copies every
        attribute before saving it while False tries to save every parameter
        as reference. Note however that it is not always possible to save the
        input as reference.
        Default is False.
    """

    plotter = PlotterDescriptor(default_type=SJIPlotter)

    def __init__(
        self,
        data,
        wcs,
        *,
        uncertainty=None,
        unit=None,
        meta=None,
        mask=None,
        copy=False,
        **kwargs,
    ) -> None:
        self.dust_masked = False
        super().__init__(
            data,
            wcs,
            uncertainty=uncertainty,
            mask=mask,
            meta=meta,
            unit=unit,
            copy=copy,
            **kwargs,
        )

    def __repr__(self) -> str:
        return f"{object.__repr__(self)}\n{self!s}"

    def __str__(self) -> str:
        if self.wcs.world_n_dim == 2:
            instance_start = self.global_coords.get("Time (UTC)")
            instance_end = None
        else:
            instance_start = self.wcs.pixel_to_world(0, 0, 0)[-1]
            instance_end = self.wcs.pixel_to_world(0, 0, self.data.shape[0] - 1)[-1]
        return textwrap.dedent(
            f"""
            SJICube
            -------
            Observatory:           {self.meta.get("TELESCOP", "IRIS")}
            Instrument:            {self.meta.get("INSTRUME")}
            Bandpass:              {self.meta.get("TWAVE1")}
            Obs Date:              {instance_start} -- {instance_end}
            Total Frames in Obs:   {self.meta.get("NBFRAMES")}
            Obs ID:                {self.meta.get("OBSID")}
            Obs Description:       {self.meta.get("OBS_DESC")}
            Axis Types:            {self.array_axis_physical_types}
            Roll:                  {self.meta.get("SAT_ROT")}
            Cube dimensions:       {self.shape}
            """,
        )

    def apply_dust_mask(self, *, undo=False):
        """
        Applies or undoes an update of the mask with the dust particles positions.

        Rewrite self.mask with/without the dust positions.

        Parameters
        ----------
        undo: `bool`
            If False, dust particles positions mask will be applied.
            If True, dust particles positions mask will be removed.
            Default=False
        """
        check_scaled(self)
        dust_mask = calculate_dust_mask(self.data)
        if undo:
            if self.mask is not None:
                self.mask = self.mask & ~dust_mask
            self.dust_masked = False
        else:
            if self.mask is None:
                self.mask = np.zeros(self.shape, dtype=bool)
            self.mask = self.mask | dust_mask
            self.dust_masked = True

    remove_cosmic_rays = remove_cosmic_rays
    remove_dust = _remove_dust

    def _get_axis_coord_index(self, axis_name, coord_loc):
        # The exposure time of a slit-jaw cube is an extra coordinate, whose axis sunraster
        # gives as a pixel axis int, where apply_exposure_time_correction unpacks a tuple
        # of array axes, as for a coordinate in the meta.
        index = super()._get_axis_coord_index(axis_name, coord_loc)
        if coord_loc == "extra_coords":
            return (self.data.ndim - 1 - index,)
        return index

    celestial_frame = property(_celestial_frame_from_cube)

    @cached_property
    def _file_frames(self):
        """
        The file's FITS header of each frame of this cube (one header for a single
        frame, a list for a time axis) and the cube's spatial item in numpy order;
        `None` when the cube has no frame headers.
        """
        headers = self.meta.get("frame_wcs_headers") if self.meta is not None else None
        low = self.wcs.low_level_wcs
        item = [slice(None)] * low.pixel_n_dim
        if isinstance(low, SlicedLowLevelWCS):
            low, item = low._wcs, list(low._slices_array)
        # Any other low-level WCS (e.g. a rebinned one) no longer matches the frame headers.
        if headers is None or not isinstance(low, gwcs.WCS):
            return None
        return headers[item[0]], tuple(item[1:])

    @cached_property
    def _frame_headers(self):
        """
        The FITS header of each frame of this cube, following its slicing: one header
        for a single frame, a list for a time axis, `None` when the cube has no frame
        headers or an integer spatial index left a 1-D image, which no header describes.
        """
        if self._file_frames is None or not all(isinstance(axis, slice) for axis in self._file_frames[1]):
            return None
        frames, spatial = self._file_frames
        if isinstance(frames, MetaDict):
            return _slice_header(frames, spatial)
        return [_slice_header(h, spatial) for h in frames]

    @cached_property
    def fits_wcs(self):
        """
        Returns a standard WCS instead of gWCS.

        One per frame, following the cube's slicing; after an integer spatial index it
        is the sliced WCS of each frame.
        """
        if self._file_frames is None:
            return None
        frames, spatial = self._file_frames
        if isinstance(frames, MetaDict):
            return WCS(frames).slice(spatial, numpy_order=True)
        return [WCS(h).slice(spatial, numpy_order=True) for h in frames]

    def to_maps(self, index: int | list[int] | None = None):
        """
        Return SunPy Maps for the requested frame(s).

        Parameters
        ----------
        index : int, list, optional
            The index of the SJI steps you want.
            By default None which will return the entire cube as a map sequence.

        Returns
        -------
        `sunpy.map.Map` or `sunpy.map.MapSequence`
            A single Map if index is an int, otherwise a MapSequence.
            Each map is built from its frame's FITS header, so sunpy picks the map class
            from the instrument (`~sunpy.map.sources.SJIMap` for slit-jaw images) and the
            map has the data unit and the frame's own exposure time.
        """
        from sunpy.map import Map  # NOQA: PLC0415

        # The frame headers say DN, which the raw integers of a memmap read are not
        check_scaled(self)
        if self._frame_headers is None:
            msg = (
                "This cube has no FITS WCS header for its frames (for example, it was rebinned or indexed along a "
                "spatial axis), so it cannot be converted to maps."
            )
            raise ValueError(msg)
        if isinstance(index, Integral):
            idx_list = [index]
        elif index is None:
            idx_list = range(self.data.shape[0])
        else:
            idx_list = index

        is_2d = self.wcs.world_n_dim == 2
        if is_2d:
            frames = [(self.data, self._frame_headers)]
        else:
            # The data and header lists wrap negative indices; resolve them once for both.
            idx_list = [range(self.data.shape[0])[i] for i in idx_list]
            frames = [(self.data[i], self._frame_headers[i]) for i in idx_list]
        # Keep the requested order; a MapSequence sorts by date by default
        maps = Map(frames, sequence=True, sortby=None)
        # Hinode/SOT cubes and most AIA channels have no IRIS colormap
        cmap = f"irissji{int(self.meta['TWAVE1'])}"
        if cmap in mpl.colormaps:
            for m in maps:
                m.plot_settings["cmap"] = cmap
        return maps[0] if is_2d or isinstance(index, Integral) else maps


class AIACube(SJICube):
    """
    Subclass of the SJICube.

    It is the same outside of the name.
    """

    def __str__(self) -> str:
        return super().__str__().replace("SJICube", "AIACube")


class SOTCube(SJICube):
    """
    Hinode/SOT image and magnetic data supplied with IRIS observations :cite:p:`itn32`.
    """

    def __str__(self) -> str:
        return super().__str__().replace("SJICube", "SOTCube")
