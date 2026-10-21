"""
This module provides general utility functions for IRIS Responses.
"""

import warnings

import numpy as np
import scipy
import scipy.io
from erfa import ErfaWarning
from scipy.interpolate import make_interp_spline

import astropy.units as u
from astropy.time import Time
from astropy.units.quantity import Quantity

from sunpy.time import parse_time

from irispy.data import ROOTDIR

__all__ = ["get_interpolated_effective_area", "get_latest_response"]


def get_latest_response(
    observation_time=None,
):
    """
    Returns the latest IRIS response structure.

    This is not downloading the latest response file from the IRIS website,'
    but rather reading the latest response file from the local data directory.

    Goal is to replicate the base functionality of the IDL routine
    ``iris_get_response.pro`` in the SSWIDL package.

    There are no plans to support anything but the latest response file.

    This routine does calculate time dependent effective areas using
    as is done in the SSWIDL version of this code.

    As a result of being translated from IDL, it has that IDL code smell and layout.
    The upside to this is that if upstream IDL changes occur, it is easy to adapt this
    in the future compared to if it was Python-ized.

    Parameters
    ----------
    observation_time : `astropy.time.Time`, `str` or `datetime.datetime`, optional
        The observation time, in any format `sunpy.time.parse_time` accepts.
        One time only; call this function once for each time.
        If not provided, the current time is used.

    Returns
    -------
    `dict`
        The keys are the names of the IDL response structure, in upper case:

        ``DATE_OBS`` : `str`
            The observation time, in ISO format.
        ``LAMBDA`` : `astropy.units.Quantity`
            Wavelength grid of the effective areas, shape ``(3601,)``.
        ``AREA_SG``, ``AREA_SJI`` : `astropy.units.Quantity`
            Effective areas in square centimeters at the observation time, shape ``(2, 3601)`` for the
            spectrograph (FUV, NUV) and ``(4, 3601)`` for the slit-jaw channels.
        ``NAME_SG``, ``NAME_SJI`` : `numpy.ndarray`
            Names of the spectrograph and slit-jaw channels.
        ``DN2PHOT_SG``, ``DN2PHOT_SJI`` : `numpy.ndarray`
            Photons per DN for each of those channels.
        ``COEFFS_FUV``, ``COEFFS_NUV``, ``COEFFS_SJI`` : `numpy.ndarray`
            Coefficients of the best-fit throughput function of time, one set per interval of
            ``C_F_TIME``, ``C_N_TIME`` and ``C_S_TIME`` (`astropy.time.Time`) and, for the spectrograph, per
            wavelength ``C_F_LAMBDA`` and ``C_N_LAMBDA`` (`astropy.units.Quantity`).
        ``GEOM_AREA`` : `astropy.units.Quantity`
            Geometric collecting area.
        ``ELEMENTS``, ``INDEX_EL_SG``, ``INDEX_EL_SJI`` : `numpy.recarray`, `numpy.ndarray`
            The optical elements and which of them each channel uses.
        ``COMMENT`` : `bytes`
            Comment of the response file.
        ``VERSION`` : `int`
            Version of the response file.
        ``VERSION_DATE`` : `astropy.time.Time`
            Date of that version.

    Raises
    ------
    ValueError
        If ``observation_time`` is not a single time, such as an array of one or more times.

    Examples
    --------
    >>> from irispy.utils.response import get_latest_response
    >>> response = get_latest_response("2020-01-01")
    >>> response["AREA_SG"].shape
    (2, 3601)
    >>> response["VERSION"]
    9
    """
    observation_time = parse_time("now" if observation_time is None else observation_time)
    if not observation_time.isscalar:
        msg = f"observation_time must be a single time, not an array of shape {observation_time.shape}"
        raise ValueError(msg)
    raw_response_data = scipy.io.readsav(ROOTDIR / "iris_sra_c_20231106.geny", python_dict=True)
    iris_response = _record_to_dict(raw_response_data["p0"])
    iris_response = {
        k: val[0] if isinstance(val, np.ndarray) and val.ndim == 1 else val for k, val in iris_response.items()
    }
    # Convert to appropriate types
    iris_response["VERSION"] = int(iris_response["VERSION"])
    iris_response["LAMBDA"] = Quantity(iris_response["LAMBDA"], unit=u.nm)
    iris_response["GEOM_AREA"] = Quantity(iris_response["GEOM_AREA"], unit=u.cm**2)
    iris_response["C_F_TIME"] = Time(iris_response["C_F_TIME"], format="utime")
    iris_response["C_F_LAMBDA"] = Quantity(iris_response["C_F_LAMBDA"], unit="nm")
    iris_response["C_N_TIME"] = Time(iris_response["C_N_TIME"], format="utime")
    iris_response["C_N_LAMBDA"] = Quantity(iris_response["C_N_LAMBDA"], unit="nm")
    iris_response["C_S_TIME"] = Time(iris_response["C_S_TIME"], format="utime")
    iris_response["ELEMENTS"]["DATE"] = parse_time(iris_response["ELEMENTS"]["DATE"].astype(str).tolist())
    iris_response["VERSION_DATE"] = parse_time(iris_response["VERSION_DATE"])
    # The effective areas at the observation time, in cm2, are filled in below
    iris_response["AREA_SG"] = np.zeros(iris_response["AREA_SG"].shape)
    iris_response["AREA_SJI"] = np.zeros(iris_response["AREA_SJI"].shape)
    iris_response["DATE_OBS"] = observation_time.isot

    # FUV SG effective areas
    lambda_range_fuv = np.array([[133.1, 135.9], [138.8, 140.8]])
    # Rough SG spectral ranges.
    # Setting effective area to 0 outside of these.
    # Time-dependent response at each of the 3 wavelengths of C_F_LAMBDA
    iris_fit_fuv = np.array(
        [_fit_xput_lite(observation_time, iris_response["C_F_TIME"], coeffs) for coeffs in iris_response["COEFFS_FUV"]]
    )
    # Interpolate onto lambda grid, separately for each of the two FUV CCD's.
    for j in range(2):
        w_fuv = np.logical_and(
            iris_response["LAMBDA"].value >= lambda_range_fuv[j, 0],
            iris_response["LAMBDA"].value <= lambda_range_fuv[j, 1],
        )
        # For some reason in Python the above for Si IV misses the last point
        # so I add an extra True value in the mask
        if j == 1:
            idx = np.argwhere(w_fuv)[-1] + 1
            w_fuv[idx] = True
        iris_response["AREA_SG"][0, w_fuv] = make_interp_spline(
            iris_response["C_F_LAMBDA"][j : j + 2],
            iris_fit_fuv[j : j + 2],
            k=1,
        )(iris_response["LAMBDA"][w_fuv])

    # The IDL code has a special case here for version 9 and later.
    # Remove wavelength dependence for the Si IV part of the FUV window.
    w_fuv_si_iv = np.logical_and(
        iris_response["LAMBDA"].value >= lambda_range_fuv[1, 0],
        iris_response["LAMBDA"].value <= lambda_range_fuv[1, 1],
    )
    # For some reason in Python the above for Si IV misses the last point
    # so I add an extra True value in the mask
    idx = np.argwhere(w_fuv_si_iv)[-1] + 1
    w_fuv_si_iv[idx] = True
    iris_response["AREA_SG"][0, w_fuv_si_iv] = np.mean(iris_response["AREA_SG"][0, w_fuv_si_iv])

    # NUV SG effective areas
    lambda_range_nuv = np.array([278.2, 283.5])
    # Rough SG spectral ranges.
    # Setting effective area to 0 outside of these.
    # Time-dependent response at each wavelength of C_N_LAMBDA
    iris_fit_nuv = np.array(
        [_fit_xput_lite(observation_time, iris_response["C_N_TIME"], coeffs) for coeffs in iris_response["COEFFS_NUV"]]
    )
    # Interpolate onto lambda grid
    w_nuv = np.where(
        np.bitwise_and(
            iris_response["LAMBDA"].value >= lambda_range_nuv[0],
            iris_response["LAMBDA"].value <= lambda_range_nuv[1],
        ),
    )
    iris_response["AREA_SG"][1, w_nuv] = make_interp_spline(
        iris_response["C_N_LAMBDA"], iris_fit_nuv, k=3, bc_type="natural"
    )(iris_response["LAMBDA"][w_nuv])

    # SJI effective areas
    for nuv in range(2):
        # Calculate baseline SJI area curves
        area_sji = iris_response["GEOM_AREA"]
        for m in range(len(iris_response["INDEX_EL_SJI"][nuv * 2])):
            index_values1 = iris_response["INDEX_EL_SJI"][nuv * 2 : nuv * 2 + 2, m]
            area_sji = area_sji * iris_response["ELEMENTS"][index_values1].trans
        # Apply time dependent profile shape adjustment to FUV SJI
        if nuv == 0:
            # FUV: apply FUV SG "slant", then normalize so that a weighted (2.4:1)
            # sum at C II and Si IV gives constant response
            weight = np.array([2.4, 1])  # Typical solar ratio CII : SiIV
            wavelength = iris_response["C_F_LAMBDA"]
            n_wavelength = len(wavelength)
            wavelength = np.array(
                [
                    wavelength[0].value,
                    ((wavelength[n_wavelength - 2] * 2 + wavelength[n_wavelength - 1]) / 3).value,
                ],
            )  # 2 wavelengths in nm
            # Calculate baseline SG area for scaling purposes
            area_sg = iris_response["GEOM_AREA"]
            for m in range(len(iris_response["INDEX_EL_SG"][nuv])):
                index_values2 = iris_response["INDEX_EL_SG"][nuv, m]
                area_sg = area_sg * iris_response["ELEMENTS"][index_values2].trans
            # SG and SJI areas at wavelength
            interpol_sg = make_interp_spline(
                iris_response["LAMBDA"],
                area_sg,
                k=1,
            )
            area_sg2 = interpol_sg(wavelength)
            area_sj2 = np.zeros((2, 2))
            for m in range(2):
                interpol_sji = make_interp_spline(
                    iris_response["LAMBDA"],
                    area_sji[m],
                    k=1,
                )
                area_sj2[:, m] = interpol_sji(wavelength)
            # Calculate the normalized slant function scal, apply to area_sji
            # Best-estimate slant, i.e., eff.area @ wavelength / baseline SG @ wavelength
            interpol_sg2 = make_interp_spline(
                iris_response["LAMBDA"],
                iris_response["AREA_SG"][0],
                k=1,
            )
            sca2 = interpol_sg2(wavelength) / area_sg2
            # Normalize slant so that total(wei*asj2*sca2)/total(wei*asj2)=1
            for m in range(2):
                sca2n = sca2 * np.sum(weight * area_sj2[:, m]) / np.sum(weight * area_sj2[:, m] * sca2)
                interpol_sca = make_interp_spline(
                    wavelength,
                    sca2n,
                    k=1,
                )
                sca1n = interpol_sca(iris_response["LAMBDA"])
                sca1n = np.clip(sca1n, a_min=0, a_max=None)
                iris_response["AREA_SJI"][m] = area_sji[m] * sca1n
        else:
            # NUV is not processed
            iris_response["AREA_SJI"][2:4] = list(area_sji)
    for j in range(4):
        # SJI specific time dependency
        iris_response["AREA_SJI"][j] *= _fit_xput_lite(
            observation_time,
            iris_response["C_S_TIME"][j, :],
            iris_response["COEFFS_SJI"][j, :],
        )
    iris_response["AREA_SG"] = Quantity(iris_response["AREA_SG"], unit=u.cm**2)
    iris_response["AREA_SJI"] = Quantity(iris_response["AREA_SJI"], unit=u.cm**2)
    return iris_response


def _fit_xput_lite(observation_time, time_cal_coeffs, cal_coeffs):
    """
    To calculate the coefficients of best-fit time function for throughput, for which we
    apply a fit based on ``cal_coeffs``.

    The procedure involved in this function is as follows:

    1. The time difference (in years) is computed from the ``observation_time`` and ``time_cal_coeffs``.
    2. A least-squares fit is performed to determine the best fit for the time-dependent effective areas given the time difference.

    Goal is to replicate the base functionality of the IDL routine
    ``fit_iris_xput.pro`` in the SSWIDL package but without the optional keyword argument.

    Parameters
    ----------
    observation_time: `astropy.time.Time`
        Observation time, a single time.
    time_cal_coeffs: `astropy.time.Time`
        Start and end times of intervals of constant ``cal_coeffs[i]``.
        These should be in "utime" format.
    cal_coeffs: `numpy.ndarray`
        Coefficients of best-fit function, with at least two columns.

    Returns
    -------
    `float`
        The fit used to compute the effective area at ``observation_time``.
    """
    if time_cal_coeffs.shape[1] != 2 or cal_coeffs.shape[1] < 2:
        # Raise ValueError as time coefficient have the wrong format.
        msg = "Incorrect number of elements either in time_cal_coeffs or in cal_coeffs."
        raise ValueError(msg)
    # Some time transformations.
    # Exponent for transition between exp.decay intervals.
    transition_exp = 1.5
    aux_cal_coeffs = np.zeros(2 * time_cal_coeffs.shape[0])
    # Looking for the closest time in the calibration time intervals.
    # Differences are given in years before passing to the next stage.
    # The last interval ends in 2049, past ERFA's leap-second table, which gives a "dubious year" warning.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ErfaWarning)
        t_diff = observation_time - time_cal_coeffs
    t_diff = t_diff.flatten().to_value(u.year)
    idx = np.where(t_diff < 0)[0]
    if idx.size == 0:
        idx = 1
    else:
        idx = idx[0]
        if idx == 0:
            idx = 1
    # If the observation_time is between the calibration time intervals of a
    # calibration file (idx % !=0) then the aux_coeffs are given by an
    # exponential (coefficient and exponential value).
    # If the observation_time is between the end calibration time interval of
    # a calibration file (cal_file_t) and the beginning calibration
    # time interval of the next calibration file (cal_file_t+1)
    # (idx% 2 == 0) then, the aux_coeffs are given by 4 values
    # corresponding to a partial exponential obtained from
    # cal_file_t and a complementary exponential obtained from the
    # cal_file_t+1
    if idx % 2 != 0:  # I.e., if idx is not even...
        dtt_0 = 1
        exp_0 = np.exp(cal_coeffs[idx // 2, 2] * (t_diff[idx - 1]))
        aux_cal_coeffs[idx - 1 : idx + 1] = np.array([dtt_0, dtt_0 * exp_0])
    else:
        dtt_1 = (t_diff[idx - 1] / (t_diff[idx - 1] - t_diff[idx])) ** transition_exp
        dtt_0 = 1 - dtt_1
        exp_0 = np.exp(cal_coeffs[(idx // 2) - 1, 2] * (t_diff[idx - 2]))
        exp_1 = np.exp(cal_coeffs[idx // 2, 2] * (t_diff[idx]))
        aux_cal_coeffs[idx - 2 : idx + 2] = np.array([dtt_0, dtt_0 * exp_0, dtt_1, dtt_1 * exp_1])
    return np.matmul(aux_cal_coeffs, cal_coeffs[:, :2].reshape(aux_cal_coeffs.shape[0]))


def get_interpolated_effective_area(iris_response, detector_type, obs_wavelength):
    """
    To compute the interpolated time-dependent effective area.

    It will generalize to the time of the observation.

    Parameters
    ----------
    iris_response : dict
        The IRIS response data loaded from `irispy.utils.response.get_latest_response`.
    detector_type : `str`
        Detector type: 'FUV' or 'NUV'.
    obs_wavelength : `astropy.units.Quantity`
        The wavelength at which the observation has been taken in Angstroms.

    Returns
    -------
    `astropy.units.Quantity`
        The effective area(s) determined by interpolation with a spline fit.
        Wavelengths outside the spectral ranges covered by the response file
        are returned as NaN.

    Notes
    -----
    The response file only defines the effective area within the nominal
    spectral ranges of each band (e.g., ~1388-1408 Å for FUV2); outside of
    these, the upstream IDL code (``iris_get_response.pro``) sets it to zero.
    Some observations read out spectral windows wider than these nominal
    ranges, and dividing by a zero (or partially interpolated near-zero)
    effective area there produces infinite or wildly amplified values, as the
    IDL code (``iris_calib_spectrum.pro``) silently does. Instead, this
    function returns NaN for such wavelengths, making it explicit in the data
    that no calibration is available there.
    """
    if detector_type.startswith("FUV"):
        detector_type_index = 0
    elif detector_type.startswith("NUV"):
        detector_type_index = 1
    else:
        msg = "Detector type not recognized."
        raise ValueError(msg)
    eff_area = iris_response["AREA_SG"][detector_type_index, :]
    response_wavelength = iris_response["LAMBDA"]
    obs_wavelength = obs_wavelength.to(u.nm)
    # Interpolate the effective areas to cover the wavelengths at which the data is recorded
    tck = make_interp_spline(
        response_wavelength,
        eff_area,
        k=1,
    )
    eff_area_interp = tck(obs_wavelength)
    # The response file sets the effective area to 0 outside the nominal spectral
    # ranges. Only wavelengths bracketed by two nonzero response points are
    # covered; interpolating the support indicator keeps the gap between the two
    # FUV CCDs excluded as well.
    supported = make_interp_spline(response_wavelength, (eff_area > 0).astype(float), k=1)(obs_wavelength) > 0.999
    eff_area_interp[~supported] = np.nan
    return eff_area_interp * u.cm**2


def _record_to_dict(arr: np.ndarray):
    """
    Convert a structured numpy array to a regular dictionary.
    """
    new_arr = arr.copy()
    for name in arr.dtype.names:
        if isinstance(arr[name], np.ndarray) and arr[name].dtype == np.object_:
            new_arr[name] = arr[name].tolist()
    return {name: new_arr[name] for name in new_arr.dtype.names}
