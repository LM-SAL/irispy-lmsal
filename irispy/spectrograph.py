import textwrap

import astropy.units as u
import numpy as np
from astropy.time import Time

from sunraster import SpectrogramCube, SpectrogramSequence
from sunraster.spectrogram import APPLY_EXPOSURE_TIME_ERROR

from irispy import utils

__all__ = ["IRISSpectrogramCube", "IRISSpectrogramCubeSequence"]


class IRISSpectrogramCube(SpectrogramCube):
    """
    Class representing IRISSpectrogramCube data described by a single WCS.

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
        mask=None,
        copy=False,
    ):
        required_meta_keys = ["detector type"]
        if not all([key in list(meta) for key in required_meta_keys]):
            raise ValueError(f"Meta must contain following keys: {required_meta_keys}")
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
        return IRISSpectrogramCube(
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
        roll = self.meta.get("SAT_ROT", None)
        if self.global_coords and "time" in self.global_coords:
            instance_start = self.global_coords["time"].min().isot
            instance_end = self.global_coords["time"].max().isot
        elif self.extra_coords and isinstance(self.axis_world_coords("time", wcs=self.extra_coords)[0], Time):
            instance_start = self.axis_world_coords("time", wcs=self.extra_coords)[0].min().isot
            instance_end = self.axis_world_coords("time", wcs=self.extra_coords)[0].max().isot
        else:
            instance_start = None
            instance_end = None
        return textwrap.dedent(
            f"""
                IRISSpectrogramCube
                -------------------
                {utils.produce_obs_repr_string(self.meta)}

                Spectrogram period: {instance_start} -- {instance_end}
                Data shape: {self.dimensions}
                Axis Types: {self.array_axis_physical_types}
                Roll: {roll}
                """
        )

    def convert_to(self, new_unit_type, time_obs=None, response_version=4):
        """
        Converts data, unit and uncertainty attributes to new unit type.

        Takes into consideration also the observation time and response version.

        The presence or absence of the exposure time correction is
        preserved in the conversions.

        Parameters
        ----------
        new_unit_type: `str`
           Unit type to convert data to. Three values are accepted:
           "DN": Relevant IRIS data number based on detector type.
           "photons": photon counts
           "radiance": Perorms radiometric calibration conversion.
        time_obs: an `astropy.time.Time` object, as a kwarg, valid for version > 2, optional
           Observation times of the datapoints.
           Must be in the format of, e.g.,
           time_obs=Time('2013-09-03', format='utime'),
           which yields 1094169600.0 seconds in value.
           The argument time_obs is ignored for versions 1 and 2.
        response_version: `int`, optional
            Version number of effective area file to be used. Cannot be set
            simultaneously with response_file or pre_launch kwarg. Default=4.

        Returns
        -------
        `IRISSpectrogramCube`
            New IRISSpectrogramCube in new units.
        """
        detector_type = utils.get_detector_type(self.meta)
        time_obs = time_obs
        response_version = response_version  # Should default to latest
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
            # TODO: spectral_data_index UNUSED
            spectral_data_index = (-1) * (np.arange(len(self.dimensions)) + 1)[spectral_wcs_index]  # NOQA
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
                self = IRISSpectrogramCube(
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
                (self.data, self.uncertainty.array), self.unit, new_unit
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
        new_cube = IRISSpectrogramCube(
            new_data,
            self.wcs,
            new_uncertainty,
            new_unit,
            self.meta,
            mask=self.mask,
        )
        new_cube._extra_coords = self.extra_coords
        return new_cube


class IRISSpectrogramCubeSequence(SpectrogramSequence):
    """
    Class for holding, slicing and plotting IRIS spectrogram data.

    This class contains all the functionality of its super class with
    some additional functionalities.

    Parameters
    ----------
    data_list: `list`
        List of `IRISSpectrogramCube` objects from the same spectral window and OBS ID.
        Must also contain the 'detector type' in its meta attribute.
    meta: `dict` or header object, optional
        Metadata associated with the sequence.
    common_axis: `int`, optional
        The axis of the NDCubes corresponding to time.
    """

    def __init__(self, data_list, meta=None, common_axis=0):
        detector_type_key = "detector type"
        # Check that meta contains required keys.
        required_meta_keys = [
            detector_type_key,
            "spectral window",
            "brightest wavelength",
            "min wavelength",
            "max wavelength",
        ]
        if not all([key in list(meta) for key in required_meta_keys]):
            raise ValueError(f"Meta must contain following keys: {required_meta_keys}")
        # Check that all spectrograms are from same specral window and OBS ID.
        if len(np.unique([cube.meta["OBSID"] for cube in data_list])) != 1:
            raise ValueError("Constituent IRISSpectrogramCube objects must have same " "value of 'OBSID' in its meta.")
        if len(np.unique([cube.meta["spectral window"] for cube in data_list])) != 1:
            raise ValueError(
                "Constituent IRISSpectrogramCube objects must have same " "value of 'spectral window' in its meta."
            )
        # Initialize Sequence.
        super().__init__(data_list, meta=meta, common_axis=common_axis)

    def __repr__(self):
        roll = self[0].meta.get("SAT_ROT", None)
        return """
IRISSpectrogramCubeSequence
---------------------------
{obs_repr}

Sequence period: {inst_start} -- {inst_end}
Sequence Shape: {seq_shape}
Axis Types: {axis_types}
Roll: {roll}

""".format(
            obs_repr=utils.produce_obs_repr_string(self.data[0].meta),
            inst_start=self[0].extra_coords["time"]["value"][0],
            inst_end=self[-1].extra_coords["time"]["value"][-1],
            seq_shape=self.dimensions,
            axis_types=self.array_axis_physical_types,
            roll=roll,
        )

    def convert_to(self, new_unit_type, copy=False):
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
            return IRISSpectrogramCubeSequence(converted_data_list, meta=self.meta, common_axis=self._common_axis)
        else:
            self.data = converted_data_list
