"""
This module provides general utility functions called by code in spectrograph.
"""


import astropy.units as u
import numpy as np
from astropy import constants

from irispy.utils.constants import RADIANCE_UNIT
from irispy.utils.utils import get_interpolated_effective_area

__all__ = [
    "convert_between_DN_and_photons",
    "convert_or_undo_photons_per_sec_to_radiance",
    "calculate_photons_per_sec_to_radiance_factor",
    "reshape_1D_wavelength_dimensions_for_broadcast",
]


def convert_between_DN_and_photons(old_data_arrays, old_unit, new_unit):
    """
    Converts arrays from IRIS DN to photons or vice versa.

    In this function, an inverse time component due to exposure time
    correction is ignored during calculations but preserved in final unit.

    Parameters
    ----------
    old_data_arrays: iterable of `numpy.ndarray`
        Arrays of data to be converted.
    old_unit: `astropy.unit.Unit`
        Unit of data arrays.
    new_unit: `astropy.unit.Unit`
        Unit to convert data arrays to.

    Returns
    -------
    `list` of `numpy.ndarray`
        Data arrays converted to new_unit.
    `astropy.unit.Unit`
        Unit of new data arrays with any inverse time component preserved.
    """
    if old_unit == new_unit or old_unit == new_unit / u.s:
        new_data_arrays = [data for data in old_data_arrays]
        new_unit_time_accounted = old_unit
    else:
        # During calculations, the time component due to exposure
        # time correction, if it has been applied, is ignored.
        # Check here whether the time correction is present in the
        # original unit so that is carried through to new unit.
        if u.s not in (old_unit * u.s).decompose().bases:
            old_unit_without_time = old_unit * u.s
            new_unit_time_accounted = new_unit / u.s
        else:
            old_unit_without_time = old_unit
            new_unit_time_accounted = new_unit
        # Convert data and uncertainty to new unit.
        new_data_arrays = [(data * old_unit_without_time).to(new_unit).value for data in old_data_arrays]
    return new_data_arrays, new_unit_time_accounted


def convert_or_undo_photons_per_sec_to_radiance(
    data_quantities,
    time_obs,
    response_version,
    obs_wavelength,
    detector_type,
    spectral_dispersion_per_pixel,
    solid_angle,
    undo=False,
):
    """
    Converts data quantities from counts/s to radiance (or vice versa).

    Parameters
    ----------
    data_quantities: iterable of `astropy.units.Quantity`
        Quantities to be converted.  Must have units of counts/s or
        radiance equivalent counts, e.g. erg / cm**2 / s / sr / Angstrom.
    time_obs: an `astropy.time.Time` object, as a kwarg, valid for version > 2
        Observation times of the datapoints.
        Must be in the format of, e.g.,
        time_obs parse_time('2013-09-03', format='utime'),
        which yields 1094169600.0 seconds in value.
        The argument time_obs is ignored for versions 1 and 2.
    response_version : `int`
        Version number of effective area file to be used. Cannot be set
        simultaneously with response_file or pre_launch kwarg. Default=4.
    obs_wavelength: `astropy.units.Quantity`
        Wavelength at each element along spectral axis of data quantities.
    detector_type: `str`
        Detector type: 'FUV', 'NUV', or 'SJI'.
    spectral_dispersion_per_pixel: scalar `astropy.units.Quantity`
        spectral dispersion (wavelength width) of a pixel.
    solid_angle: scalar `astropy.units.Quantity`
        Solid angle corresponding to a pixel.
    undo: `bool`
        If False, converts counts/s to radiance.
        If True, converts radiance to counts/s.
        Default=False

    Returns
    -------
    `list` of `astropy.units.Quantity`
        Data quantities converted to radiance or counts/s
        depending on value of undo kwarg.
    """
    # Check data quantities are in the right units.
    if undo is True:
        for i, data in enumerate(data_quantities):
            if not data.unit.is_equivalent(RADIANCE_UNIT):
                raise ValueError(
                    "Invalid unit provided.  As kwarg undo=True, "
                    "unit must be equivalent to {}.  Error found for {}th element "
                    "of data_quantities. Unit: {}".format(RADIANCE_UNIT, i, data.unit)
                )
    else:
        for data in data_quantities:
            if data.unit != u.photon / u.s:
                raise ValueError(
                    "Invalid unit provided.  As kwarg undo=False, "
                    "unit must be equivalent to {}.  Error found for {}th element "
                    "of data_quantities. Unit: {}".format(u.photon / u.s, i, data.unit)
                )
    photons_per_sec_to_radiance_factor = calculate_photons_per_sec_to_radiance_factor(
        time_obs,
        response_version,
        obs_wavelength,
        detector_type,
        spectral_dispersion_per_pixel,
        solid_angle,
    )
    # Change shape of arrays so they are compatible for broadcasting
    # with data and uncertainty arrays.
    photons_per_sec_to_radiance_factor = reshape_1D_wavelength_dimensions_for_broadcast(
        photons_per_sec_to_radiance_factor, data_quantities[0].ndim
    )
    # Perform (or undo) radiometric conversion.
    if undo is True:
        new_data_quantities = [
            (data / photons_per_sec_to_radiance_factor).to(u.photon / u.s) for data in data_quantities
        ]
    else:
        new_data_quantities = [
            (data * photons_per_sec_to_radiance_factor).to(RADIANCE_UNIT) for data in data_quantities
        ]
    return new_data_quantities


def calculate_photons_per_sec_to_radiance_factor(
    time_obs,
    response_version,
    wavelength,
    detector_type,
    spectral_dispersion_per_pixel,
    solid_angle,
):
    """
    Calculates multiplicative factor that converts counts/s to radiance for
    given wavelengths.

    Parameters
    ----------
    time_obs: an `astropy.time.Time` object, as a kwarg, valid for version > 2
        Observation times of the datapoints.
        Must be in the format of, e.g.,
        time_obs=parse_time('2013-09-03', format='utime'),
        which yields 1094169600.0 seconds in value.
        The argument time_obs is ignored for versions 1 and 2.
    response_version : `int`
        Version number of effective area file to be used. Cannot be set
        simultaneously with response_file or pre_launch kwarg. Default=4.
    wavelength: `astropy.units.Quantity`
        Wavelengths for which counts/s-to-radiance factor is to be calculated
    detector_type: `str`
        Detector type: 'FUV' or 'NUV'.
    spectral_dispersion_per_pixel: scalar `astropy.units.Quantity`
        spectral dispersion (wavelength width) of a pixel.
    solid_angle: scalar `astropy.units.Quantity`
        Solid angle corresponding to a pixel.

    Returns
    -------
    `astropy.units.Quantity`
        Mutliplicative conversion factor from counts/s to radiance units
        for input wavelengths.
    """
    # Get effective area and interpolate to observed wavelength grid.
    eff_area_interp = get_interpolated_effective_area(
        time_obs, response_version, detector_type, obs_wavelength=wavelength
    )
    # Return radiometric conversed data assuming input data is in units of photons/s.
    return (
        constants.h
        * constants.c
        / wavelength
        / u.photon
        / spectral_dispersion_per_pixel
        / eff_area_interp
        / solid_angle
    )


def reshape_1D_wavelength_dimensions_for_broadcast(wavelength, n_data_dim):
    if n_data_dim == 1:
        pass
    elif n_data_dim == 2:
        wavelength = wavelength[np.newaxis, :]
    elif n_data_dim == 3:
        wavelength = wavelength[np.newaxis, np.newaxis, :]
    else:
        raise ValueError("IRISSpectrogram dimensions must be 2 or 3.")
    return wavelength
