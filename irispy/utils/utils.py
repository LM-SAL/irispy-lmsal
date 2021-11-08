"""
This module provides general utility functions.
"""

import astropy.units as u
import numpy as np
from astropy.modeling.models import custom_model
from scipy import interpolate, ndimage

from irispy.utils.response import get_iris_response

__all__ = [
    "get_detector_type",
    "get_interpolated_effective_area",
    "calculate_dust_mask",
    "gaussian1d_on_linear_bg",
]


@custom_model
def gaussian1d_on_linear_bg(
    x,
    amplitude=None,
    mean=None,
    standard_deviation=None,
    constant_term=None,
    linear_term=None,
):
    return amplitude * np.exp(-(((x - mean) / standard_deviation) ** 2)) + constant_term + linear_term * x


def get_detector_type(meta):
    """
    Gets the IRIS detector type from a meta dictionary.

    In this function, FUV1 and FUV2 are just assigned as FUV.

    Parameters
    ----------
    meta: dict-like
        Dictionary-like object containing entry for "detector type"

    Returns
    -------
    `str`
       Detector type.
    """
    if "FUV" in meta["detector type"]:
        detector_type = "FUV"
    else:
        detector_type = meta["detector type"]
    return detector_type


def get_interpolated_effective_area(time_obs, response_version, detector_type, obs_wavelength):
    """
    To compute the interpolated time-dependent effective area.

    Parameters
    ----------
    time_obs: an `astropy.time.Time` object, as a kwarg, valid for version > 2
        Observation times of the datapoints.
        Must be in the format of, e.g.,
        time_obs parse_time('2013-09-03', format='utime'),
        which yields 1094169600.0 seconds in value.
        The argument time_obs is ignored for versions 1 and 2.
    response_version : `int`
        Version number of effective area file to be used. Cannot be set
        simultaneously with response_file or pre_launch kwarg. Default=4.
    detector_type: `str`
        Detector type: 'FUV' or 'NUV'.
    obs_wavelength: `astropy.units.Quantity`
        The wavelength at which the observation has been taken in Angstroms.

    Returns
    -------
    `numpy.array`
        The effective area(s) determined by interpolation with a spline fit.
    """
    # Generalizing to the time of obs.
    time_obs = time_obs
    response_version = response_version
    iris_response = get_iris_response(time_obs, response_version)
    if detector_type == "FUV":
        detector_type_index = 0
    elif detector_type == "NUV":
        detector_type_index = 1
    else:
        raise ValueError("Detector type not recognized.")
    eff_area = iris_response["AREA_SG"][detector_type_index, :]
    response_wavelength = iris_response["LAMBDA"]
    # Interpolate the effective areas to cover the wavelengths
    # at which the data is recorded:
    eff_area_interp_base_unit = u.Angstrom
    tck = interpolate.splrep(
        response_wavelength.to(eff_area_interp_base_unit).value,
        eff_area.to(eff_area_interp_base_unit ** 2).value,
        s=0,
    )
    eff_area_interp = interpolate.splev(obs_wavelength.to(eff_area_interp_base_unit).value, tck) * (
        eff_area_interp_base_unit ** 2
    )
    return eff_area_interp


def calculate_dust_mask(data_array):
    """
    Calculate a mask with the dust positions in a given arrayself.

    Parameters
    ----------
    data_array : `numpy.ndarray`
        This array contains some dust poisition that will be calculated. The array
        must have scaled values.

    Returns
    -------
    `numpy.ndarray` of `bool`
        This array has the same shape than data_array and contains the dust positions
        when the value is True.
    """
    # Creating a mask with the same shape than the inputed data array.
    mask = np.zeros_like(data_array, dtype=bool)
    # Set the pixel value to True is the pixel is recognized as a dust pixel.
    mask[(data_array < 0.5) & (data_array > -200)] = True
    # Extending the mask to avoid the neighbours pixel influenced by the dust pixels.
    struct = np.array([np.zeros((3, 3)), np.ones((3, 3)), np.zeros((3, 3))], dtype=bool)
    mask = ndimage.binary_dilation(mask, structure=struct).astype(mask.dtype)
    return mask
