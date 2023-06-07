import logging
import textwrap

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.time import Time
from ndcube import NDCollection
from sunpy.coordinates import Helioprojective
from sunraster import SpectrogramCube as SpecCube
from sunraster import SpectrogramSequence as SpecSeq
from sunraster.meta import Meta, SlitSpectrographMetaABC
from sunraster.spectrogram import APPLY_EXPOSURE_TIME_ERROR

from irispy import utils
from irispy.utils.constants import SPECTRAL_BAND
from irispy.visualization import Plotter, _set_axis_colors

__all__ = ["Collection", "SpectrogramCube", "SpectrogramCubeSequence", "SGMeta"]


class Collection(NDCollection):
    def __str__(self):
        return textwrap.dedent(
            f"""
            Collection
            ----------
            Cube keys: {tuple(self.keys())}
            Number of Cubes: {len(self)}
            Aligned dimensions: {self.aligned_dimensions}
            Aligned physical types: {self.aligned_axis_physical_types}
            """,
        )


class SpectrogramCube(SpecCube):
    """
    Class representing SpectrogramCube data described by a single WCS.

    Parameters
    ----------
    data: `numpy.ndarray`
        The array holding the actual data in this object.
    wcs: `ndcube.wcs.wcs.WCS`
        The WCS object containing the axes' information
    unit : `astropy.unit.Unit` or `str`
        Unit for the dataset. Strings that can be converted to a Unit are allowed.
    meta : dict-like object
        Additional meta information about the dataset. Must contain at least the
        following keys:
        - detector type: str, (FUV1, FUV2 or NUV)
        - OBSID: int
        - spectral window: str
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
    """

    def __init__(
        self,
        data,
        wcs,
        uncertainty,
        unit,
        meta,
        *,
        mask=None,
        copy=False,
    ):
        super().__init__(
            data,
            wcs,
            unit=unit,
            uncertainty=uncertainty,
            mask=mask,
            meta=meta,
            copy=copy,
        )

    def __getitem__(self, item):
        result = super().__getitem__(item)
        return SpectrogramCube(
            result.data,
            result.wcs,
            result.uncertainty,
            result.unit,
            result.meta,
            mask=result.mask,
        )

    def __repr__(self):
        return f"{object.__repr__(self)}\n{str(self)}"

    def __str__(self):
        instance_start = None
        instance_end = None
        if self.global_coords and "time" in self.global_coords:
            instance_start = self.global_coords["time"].min().isot
            instance_end = self.global_coords["time"].max().isot
        elif self.extra_coords and self.axis_world_coords("time", wcs=self.extra_coords):
            instance_start = self.axis_world_coords("time", wcs=self.extra_coords)[0].min().isot
            instance_end = self.axis_world_coords("time", wcs=self.extra_coords)[0].max().isot
        return textwrap.dedent(
            f"""
            SpectrogramCube
            ---------------
            OBS ID:             {self.meta.get("OBSID")}
            OBS Description:    {self.meta.get("OBS_DESC")}
            OBS period:         {self.meta.get("STARTOBS")} -- {self.meta.get("ENDOBS")}
            Spectrogram period: {instance_start} -- {instance_end}
            Data shape:         {self.dimensions}
            Axis Types:         {self.array_axis_physical_types}
            Roll:               {self.meta.get("SAT_ROT")}
            """,
        )

    def plot(self, *args, **kwargs):
        cmap = kwargs.get("cmap")
        if not cmap:
            try:
                cmap = plt.get_cmap(name=f"irissji{int(self.meta.detector[:3])}")
            except Exception as e:  # NOQA: BLE001
                logging.debug(e)
                cmap = "viridis"
        kwargs["cmap"] = cmap
        if len(self.dimensions) == 1:
            kwargs.pop("cmap")
        ax = Plotter(ndcube=self).plot(*args, **kwargs)
        _set_axis_colors(ax)
        return ax

    def convert_to(self, new_unit_type, time_obs=None, response_version=4):
        """
        Converts data, unit and uncertainty attributes to new unit type.

        Takes into consideration also the observation time and response version.

        The presence or absence of the exposure time correction is
        preserved in the conversions.

        Parameters
        ----------
        new_unit_type : `str`
            Unit type to convert data to. Three values are accepted:
            "DN": Relevant IRIS data number based on detector type.
            "photons": photon counts
            "radiance": Perorms radiometric calibration conversion.
        time_obs : `astropy.time.Time`, optional
            Observation times of the datapoints.
            Must be in the format of, e.g.,
            ``time_obs=Time('2013-09-03')``,
            The argument time_obs is ignored for versions 1 and 2.
        response_version : `int`, optional
            Version number of effective area file to be used, by default = 6.

        Returns
        -------
        `SpectrogramCube`
            New SpectrogramCube in new units.
        """
        detector_type = utils.get_detector_type(self.meta)
        if new_unit_type == "radiance" or self.unit.is_equivalent(utils.RADIANCE_UNIT):
            # Get spectral dispersion per pixel.
            spectral_wcs_index = np.where(np.array(self.wcs.wcs.ctype) == "WAVE")[0][0]
            spectral_dispersion_per_pixel = (
                self.wcs.wcs.cdelt[spectral_wcs_index] * self.wcs.wcs.cunit[spectral_wcs_index]
            )
            # Get solid angle from slit width for a pixel.
            lat_wcs_index = ["HPLT" in c for c in self.wcs.wcs.ctype]
            lat_wcs_index = np.arange(len(self.wcs.wcs.ctype))[lat_wcs_index]
            lat_wcs_index = lat_wcs_index[0]
            solid_angle = self.wcs.wcs.cdelt[lat_wcs_index] * self.wcs.wcs.cunit[lat_wcs_index] * utils.SLIT_WIDTH
            # Get wavelength for each pixel.
            obs_wavelength = self.axis_world_coords(2)

        if new_unit_type == "DN" or new_unit_type == "photons":
            if self.unit.is_equivalent(utils.RADIANCE_UNIT):
                # Convert from radiance to counts/s
                new_data_quantities = utils.convert_or_undo_photons_per_sec_to_radiance(
                    (self.data * self.unit, self.uncertainty.array * self.unit),
                    time_obs,
                    response_version,
                    obs_wavelength,
                    detector_type,
                    spectral_dispersion_per_pixel,
                    solid_angle,
                    undo=True,
                )
                new_data = new_data_quantities[0].value
                new_uncertainty = new_data_quantities[1].value
                new_unit = new_data_quantities[0].unit
                self = SpectrogramCube(
                    new_data,
                    self.wcs,
                    new_uncertainty,
                    new_unit,
                    self.meta,
                    mask=self.mask,
                )
                self._extra_coords = self.extra_coords
            if new_unit_type == "DN":
                new_unit = utils.DN_UNIT[detector_type]
            else:
                new_unit = u.photon
            new_data_arrays, new_unit = utils.convert_between_DN_and_photons(
                (self.data, self.uncertainty.array),
                self.unit,
                new_unit,
            )
            new_data = new_data_arrays[0]
            new_uncertainty = new_data_arrays[1]
        elif new_unit_type == "radiance":
            if self.unit.is_equivalent(utils.RADIANCE_UNIT):
                new_data = self.data
                new_uncertainty = self.uncertainty
                new_unit = self.unit
            else:
                # Ensure spectrogram is in units of counts/s.
                cube = self.convert_to("photons")
                try:
                    cube = cube.apply_exposure_time_correction()
                except ValueError(APPLY_EXPOSURE_TIME_ERROR):
                    pass
                # Convert to radiance units.
                new_data_quantities = utils.convert_or_undo_photons_per_sec_to_radiance(
                    (cube.data * cube.unit, cube.uncertainty.array * cube.unit),
                    time_obs,
                    response_version,
                    obs_wavelength,
                    detector_type,
                    spectral_dispersion_per_pixel,
                    solid_angle,
                )
                new_data = new_data_quantities[0].value
                new_uncertainty = new_data_quantities[1].value
                new_unit = new_data_quantities[0].unit
        else:
            raise ValueError("Input unit type not recognized.")
        new_cube = SpectrogramCube(
            new_data,
            self.wcs,
            new_uncertainty,
            new_unit,
            self.meta,
            mask=self.mask,
        )
        new_cube._extra_coords = self.extra_coords
        return new_cube


class SpectrogramCubeSequence(SpecSeq):
    """
    Class for holding, slicing and plotting IRIS spectrogram data.

    This class contains all the functionality of its super class with
    some additional functionalities.

    Parameters
    ----------
    data_list: `list`
        List of `SpectrogramCube` objects from the same spectral window and OBS ID.
    meta: `dict` or header object, optional
        Metadata associated with the sequence.
    common_axis: `int`, optional
        The axis of the NDCubes corresponding to time.
    """

    def __init__(self, data_list, meta=None, common_axis=0):
        # Check that all spectrograms are from same spectral window and OBS ID.
        if len(np.unique([cube.meta["OBSID"] for cube in data_list])) != 1:
            raise ValueError("Constituent SpectrogramCube objects must have same value of 'OBSID' in its meta.")
        super().__init__(data_list, meta=meta, common_axis=common_axis)

    def __str__(self):
        # Overload it get the class name in the string
        return super().__str__()

    def convert_to(self, new_unit_type, *, copy=False):
        """
        Converts data, uncertainty and unit of each spectrogram in sequence to
        new unit.

        Parameters
        ----------
        new_unit_type: `str`
           Unit type to convert data to.  Three values are accepted:
           "DN": Relevant IRIS data number based on detector type.
           "photons": photon counts
           "radiance": Perorms radiometric calibration conversion.
        copy: `bool`
            If True a new instance with the converted data values is return.
            If False, the current instance is overwritten.
            Default=False
        """
        converted_data_list = []
        for cube in self.data:
            converted_data_list.append(cube.convert_to(new_unit_type))
        if copy is True:
            return SpectrogramCubeSequence(converted_data_list, meta=self.meta, common_axis=self._common_axis)
        self.data = converted_data_list
        return None


class SGMeta(Meta, metaclass=SlitSpectrographMetaABC):
    def __init__(self, header, spectral_window, **kwargs):
        super().__init__(header, **kwargs)
        spectral_windows = np.array([self[f"TDESC{i}"] for i in range(1, self["NWIN"] + 1)])
        window_mask = np.array([spectral_window in window for window in spectral_windows])
        if window_mask.sum() < 1:
            raise ValueError(
                "Spectral window not found. "
                f"Input spectral window: {spectral_window}; "
                f"Spectral windows in header: {spectral_windows}",
            )
        if window_mask.sum() > 1:
            raise ValueError(
                "Spectral window must be unique. "
                f"Input spectral window: {spectral_window}; "
                f"Ambiguous spectral windows in header: {spectral_windows[window_mask]}",
            )
        self._iwin = np.arange(len(spectral_windows))[window_mask][0] + 1

    def __str__(self):
        return textwrap.dedent(
            f"""
                SGMeta
                ------
                Observatory:     {self.observatory}
                Instrument:      {self.instrument}
                Detector:        {self.detector}
                Spectral Window: {self.spectral_window}
                Spectral Range:  {self.spectral_range}
                Spectral Band:   {self.spectral_band}
                Dimensions:      {self.dimensions}
                Date:            {self.date_reference}
                OBS ID:          {self.observing_mode_id}
                OBS Description: {self.observing_mode_description}
                """,
        )

    def __repr__(self):
        return f"{object.__repr__(self)}\n{str(self)}"

    def _construct_time(self, key):
        val = self.get(key)
        if val is not None:
            val = Time(val, format="fits", scale="utc")
        return val

    @property
    def spectral_window(self):
        return self.get(f"TDESC{self._iwin}")

    @property
    def detector(self):
        return self.get(f"TDET{self._iwin}")

    @property
    def instrument(self):
        return self.get("INSTRUME")

    @property
    def observatory(self):
        return self.get("TELESCOP")

    @property
    def processing_level(self):
        return self.get("DATA_LEV")

    @property
    def distance_to_sun(self):
        return self.get("DSUN_OBS") * u.m

    @property
    def date_reference(self):
        return self._construct_time("DATE_OBS")

    @property
    def date_start(self):
        return self.date_reference

    @property
    def date_end(self):
        return self._construct_time("DATE_END")

    @property
    def observing_mode_id(self):
        return int(self.get("OBSID"))

    # ---------- IRIS-specific metadata properties ----------
    @property
    def dimensions(self):
        return self.shape.tolist()

    @property
    def observing_mode_description(self):
        return self.get("OBS_DESC")

    @property
    def observing_campaign_start(self):
        """
        Start time of observing campaign.
        """
        return self._construct_time("STARTOBS")

    @property
    def observing_campaign_end(self):
        """
        End time of observing mode.
        """
        return self._construct_time("ENDOBS")

    @property
    def observation_includes_SAA(self):
        """
        Whether IRIS passed through SAA during observations.
        """
        return bool(self.get("SAA"))

    @property
    def satellite_rotation(self):
        """
        Satellite roll from solar north.
        """
        return self.get("SAT_ROT") * u.deg

    @property
    def exposure_control_triggers_in_observation(self):
        """
        Number of times automatic exposure control triggered during observing
        campaign.
        """
        return self.get("AECNOBS")

    @property
    def exposure_control_triggers_in_raster(self):
        """
        Number of times automatic exposure control was triggered during this
        raster.
        """
        return self.get("AECNRAS")

    @property
    def number_raster_positions(self):
        """
        Number of positions in raster.
        """
        self.get("NRASTERP")

    @property
    def spectral_range(self):
        """
        The spectral range of the spectral window.
        """
        return [self.get(f"TWMIN{self._iwin}"), self.get(f"TWMAX{self._iwin}")] * u.AA

    @property
    def spectral_band(self):
        """
        The spectral band of the spectral window.
        """
        return SPECTRAL_BAND[self.spectral_window]

    @property
    def raster_fov_width_y(self):
        """
        Width of the field of view of the raster in the Y (slit) direction.
        """
        return self.get("FOVY") * u.arcsec

    @property
    def raster_fov_width_x(self):
        """
        Width of the field of view of the raster in the X (rastering)
        direction.
        """
        return self.get("FOVX") * u.arcsec

    @property
    def fov_center(self):
        """
        Location of the center of the field of view.
        """
        return SkyCoord(
            Tx=self.get("XCEN"),
            Ty=self.get("YCEN"),
            unit=u.arcsec,
            frame=Helioprojective,
        )

    @property
    def automatic_exposure_control_enabled(self):
        return bool(self.get("IAECFLAG"))

    @property
    def tracking_mode_enabled(self):
        return bool(self.get("TR_MODE"))

    @property
    def observatory_at_high_latitude(self):
        """
        Whether IRIS passed through high Earth latitude during observations.
        """
        return bool(self.get("HLZ"))

    @property
    def spatial_summing_factor(self):
        """
        Number of pixels summed together in the spatial (Y/slit) direction.
        """
        return self.get("SUMSPAT")

    @property
    def spectral_summing_factor(self):
        """
        Number of pixels summed together in the spectral direction.
        """
        if "fuv" in self.detector.lower():
            return self.get("SUMSPTRF")
        return self.get("SUMSPTRN")
