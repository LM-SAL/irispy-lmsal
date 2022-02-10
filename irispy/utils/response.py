"""
This module provides general utility functions for IRIS Responses.
"""
import os.path
import datetime

import astropy.units as u
import numpy as np
import scipy
import scipy.io
from astropy.time import Time
from astropy.units.quantity import Quantity
from sunpy.time import parse_time
from sunpy.util.config import get_and_create_download_dir
from sunpy.util.net import download_file

IRIS_RESPONSE_REMOTE_PATH = "https://sohowww.nascom.nasa.gov/solarsoft/iris/response/"
RESPONSE_VERSION_FILENAMES = {
    "1": "iris_sra_20130211.geny",
    "2": "iris_sra_20130715.geny",
    "3": "iris_sra_c_20150331.geny",
    "4": "iris_sra_c_20161022.geny",
}

__all__ = ["fit_iris_xput", "get_iris_response"]


def get_iris_response(
    time_obs=None,
    pre_launch=False,
    response_file=None,
    response_version=None,
    force_download=False,
):
    """
    Returns IRIS response structure.

    One and only one of pre_launch, response_file and response_version must be set.

    Parameters
    ----------
    time_obs: an `astropy.time.Time` object, as a kwarg, valid for version > 2 ,optional
        Observation times of the datapoints.
        Must be in the format of, e.g.,
        time_obs=parse_time('2013-09-03', format='utime'),
        which yields 1094169600.0 seconds in value.
        The argument time_obs is ignored for versions 1 and 2.
    pre_launch: `bool`,optional
        Equivalent to setting response_version=2.  Cannot be set
        simultaneously with response_file kwarg. Default=False.
    response_file: `int`,optional
        Version number of effective area file to be used.  Cannot be set
        simultaneously with pre_launch kwarg.  Default=latest.
    response_version : `int`,optional
        Version number of effective area file to be used. Cannot be set
        simultaneously with response_file or pre_launch kwarg. Default=4.

    Notes
    -----
    This routine does calculate time dependent effective areas using
    version 3 and above of the response functions as is done in the SSW version
    of this code. This code has been updated to calculate time-dependent
    effective areas.

    Returns
    -------
    `dict`
        Various parameters regarding IRIS response or effective area structure.
        Includes the following keys:
        date_obs: `astropy.time.Time`
        lambda: `astropy.units.Quantity`
        area_sg: `astropy.units.Quantity`
        name_sg: `str`
        dn2phot_sg: `tuple` of length 2
        area_sji: `astropy.units.Quantity`
        name_sji: `str`
        dn2phot_sji:  `tuple` of length 4
        comment: `str`
        version: `int`
        version_date: `astropy.time.Time`
    """
    # Ensures the file exits in the path given.
    if response_file is not None:
        if not (os.path.isfile(response_file)):
            raise KeyError("Not a valid file path")

    # Ensure conflicting kwargs are not set.
    if pre_launch:
        pre_launch_set = True
    else:
        pre_launch_set = False
    if response_file:
        response_file_set = True
    else:
        response_file_set = False
    if response_version:
        response_version_set = True
    else:
        response_version_set = False
    if pre_launch_set + response_file_set + response_version_set > 1:
        raise ValueError("One and only one of kwargs pre_launch, response_file " "and response_version must be set.")
    # If pre_launch set, define response_version to 2.
    if pre_launch:
        response_version = 2
    if (response_file_set is False) and (response_version_set is False):
        response_version = 4
    # If response_file not set, define appropriate response file
    # based on version.
    if not response_file:
        try:
            response_filename = RESPONSE_VERSION_FILENAMES[str(response_version)]
        except KeyError:
            raise KeyError("Version number not recognized.")
        # Define the directory in which the response file should exist
        # to be the sunpy download directory.
        # Check response file exists in download_dir. If not, download it.
        path = download_file(
            IRIS_RESPONSE_REMOTE_PATH + response_filename,
            get_and_create_download_dir(),
            overwrite=force_download,
        )
        # Define response file as path + filename.
    # Read response file and store in a dictionary.
    raw_response_data = scipy.io.readsav(path)
    iris_response = {name: raw_response_data["p0"][name][0] for name in raw_response_data["p0"].dtype.names}
    # Convert some properties to more convenient types.
    iris_response["LAMBDA"] = Quantity(iris_response["LAMBDA"], unit=u.nm)
    iris_response["AREA_SG"] = Quantity(iris_response["AREA_SG"], unit=u.cm**2)
    iris_response["AREA_SJI"] = Quantity(iris_response["AREA_SJI"], unit=u.cm**2)
    iris_response["GEOM_AREA"] = Quantity(iris_response["GEOM_AREA"], unit=u.cm**2)
    iris_response["VERSION"] = iris_response["VERSION"]
    # Convert some properties not found in version below version 3 to
    # more convenient types.
    if int(iris_response["VERSION"]) > 2:
        # If DATE_OBS has a value, convert to `astropy.time.Time`, else set to None.
        try:
            iris_response["DATE_OBS"] = parse_time(iris_response["DATE_OBS"], format="utime")
        except Exception:
            iris_response["DATE_OBS"] = None

        time_obs = np.array([time_obs.value])

        # Convert C_F_TIME to array of time objects while conserving shape.
        c_f_time = np.empty(iris_response["C_F_TIME"].shape, dtype=object)
        for i, row in enumerate(iris_response["C_F_TIME"]):
            for j, t in enumerate(row):
                c_f_time[i][j] = parse_time(t, format="utime")
        iris_response["C_F_TIME"] = c_f_time

        # Convert C_F_LAMBDA to Quantity.
        iris_response["C_F_LAMBDA"] = Quantity(iris_response["C_F_LAMBDA"], unit="nm")

        # Convert C_N_TIME to array of time objects while
        # conserving shape.
        c_n_time = np.empty(iris_response["C_N_TIME"].shape, dtype=object)
        for i, row in enumerate(iris_response["C_N_TIME"]):
            for j, t in enumerate(row):
                c_n_time[i][j] = parse_time(t, format="utime")
        iris_response["C_N_TIME"] = c_n_time

        # Convert C_N_LAMBDA to Quantity.
        iris_response["C_N_LAMBDA"] = Quantity(iris_response["C_N_LAMBDA"], unit="nm")

        # Convert C_S_TIME to array of time objects while
        # conserving shape.
        c_s_time = np.empty(iris_response["C_S_TIME"].shape, dtype=object)
        for i, row in enumerate(iris_response["C_S_TIME"]):
            for j, column in enumerate(row):
                for k, t in enumerate(column):
                    c_s_time[i][j][k] = parse_time(t, format="utime")
        iris_response["C_S_TIME"] = c_s_time

        # Convert DATE in ELEMENTS array to array of time objects.
        for i, t in enumerate(iris_response["ELEMENTS"]["DATE"]):
            iris_response["ELEMENTS"]["DATE"][i] = parse_time(t, format="iso")
            # Convert VERSION_DATE to time object.
            iris_response["VERSION_DATE"] = parse_time(iris_response["VERSION_DATE"])

    if int(iris_response["VERSION"]) < 2:
        # Change DATE tag in data with version < 2 to VERSION_DATE to
        # be consistent with more recent versions.
        iris_response["VERSION_DATE"] = Time(
            datetime.datetime(
                int(iris_response["DATE"][0:4]),
                int(iris_response["DATE"][4:6]),
                int(iris_response["DATE"][6:8]),
            )
        )
        del iris_response["DATE"]

    if int(iris_response["VERSION"]) > 2 and time_obs is not None:
        try:
            n_time_obs = len(time_obs)
        except Exception:
            n_time_obs = 1
        iris_response["AREA_SG"] = np.zeros(iris_response["AREA_SG"].shape)
        iris_response["AREA_SJI"] = np.zeros(iris_response["AREA_SJI"].shape)

        # 1. FUV SG effective areas
        lambran_fuv = np.array([[133.1, 135.9], [138.8, 140.8]])
        # Rough SG spectral ranges.  Setting effective area to 0 outside of these.
        shp_fuv = iris_response["COEFFS_FUV"].shape
        # Time-dependent response for shp_0[0] = 3 wavelengths
        iris_fit_fuv = np.zeros((n_time_obs, shp_fuv[0]))
        for j in range(shp_fuv[0]):
            iris_fit_fuv[:, j] = fit_iris_xput(
                time_obs,
                iris_response["C_F_TIME"],
                iris_response["COEFFS_FUV"][j, :, :],
            )
        # Interpolate onto lambda grid, separately for each of the two FUV CCD's.
        for j in range(2):
            w_fuv = np.logical_and(
                iris_response["LAMBDA"].value >= lambran_fuv[j, 0],
                iris_response["LAMBDA"].value <= lambran_fuv[j, 1],
            )
            for k in range(n_time_obs):
                interpol_fuv = scipy.interpolate.interp1d(
                    iris_response["C_F_LAMBDA"][j : j + 2],
                    np.squeeze(iris_fit_fuv[k, j : j + 2]),
                    fill_value="extrapolate",
                )
                iris_response["AREA_SG"][0, w_fuv] = interpol_fuv(iris_response["LAMBDA"][w_fuv])

        # 2. NUV SG effective areas
        lambran_nuv = np.array([278.2, 283.5])
        # Rough SG spectral ranges.  Setting effective area to 0 outside of these.
        shp_nuv = iris_response["COEFFS_NUV"].shape
        # Time-dependent response for shp_1[0] wavelengths
        iris_fit_nuv = np.zeros((n_time_obs, shp_nuv[0]))
        for j in range(shp_nuv[0]):
            iris_fit_nuv[:, j] = fit_iris_xput(
                time_obs,
                iris_response["C_N_TIME"],
                iris_response["COEFFS_NUV"][j, :, :],
            )
        # Interpolate onto lambda grid
        w_nuv = np.where(
            np.logical_and(
                iris_response["LAMBDA"].value >= lambran_nuv[0],
                iris_response["LAMBDA"].value <= lambran_nuv[1],
            )
        )
        if int(iris_response["VERSION"]) <= 3:
            for k in range(n_time_obs):
                interpol_nuv = scipy.interpolate.interp1d(
                    iris_response["C_N_LAMBDA"][:],
                    np.squeeze(iris_fit_nuv[k, :]),
                    fill_value="extrapolate",
                )
                iris_response["AREA_SG"][1, w_nuv] = interpol_nuv(iris_response["LAMBDA"][w_nuv])
        else:
            for k in range(n_time_obs):
                interpol_nuv = scipy.interpolate.CubicSpline(
                    iris_response["C_N_LAMBDA"][:],
                    np.squeeze(iris_fit_nuv[k, :]),
                    extrapolate=True,
                    bc_type="natural",
                    axis=0,
                )
                iris_response["AREA_SG"][1, w_nuv] = interpol_nuv(iris_response["LAMBDA"][w_nuv])

        # 3. SJI effective areas
        if int(iris_response["VERSION"]) <= 3:  # Meaning for version 3 only
            shp_sji = iris_response["COEFFS_SJI"].shape
            for j in range(shp_sji[0]):
                # Calculate pre-launch area from the individual elements
                prelaunch_area = iris_response["GEOM_AREA"]
                for k in range(len(iris_response["INDEX_EL_SJI"][j, :])):
                    index_values0 = iris_response["INDEX_EL_SJI"][j, k]
                    prelaunch_area = prelaunch_area * iris_response["ELEMENTS"][index_values0].trans
                # Time dependent response
                iris_fit_sji = fit_iris_xput(
                    time_obs,
                    iris_response["C_S_TIME"][j, :, :],
                    iris_response["COEFFS_SJI"][j, :, :],
                )
                # Time dependetn profiles
                for k in range(n_time_obs):
                    iris_response["AREA_SJI"][j, :] = prelaunch_area * iris_fit_sji[k]
        else:  # For version 4 and above
            for nuv in range(2):
                # Calculate baseline SJI area curves
                area_sji = iris_response["GEOM_AREA"]
                for m in range(len(iris_response["INDEX_EL_SJI"][nuv * 2, :])):
                    index_values1 = iris_response["INDEX_EL_SJI"][nuv * 2 : nuv * 2 + 2, m]
                    area_sji = area_sji * iris_response["ELEMENTS"][index_values1].trans
                # Apply time dependent profile shape adjustment to FUV SJI
                if nuv == 0:
                    # FUV: apply FUV SG "slant", then normalize so that a weighted (2.4:1)
                    # sum at C II and Si IV gives constant response
                    weight = np.array([2.4, 1.0])  # Typical solar ratio CII : SiIV
                    wavelength = iris_response["C_F_LAMBDA"]
                    n_wavelength = len(wavelength)
                    wavelength = np.array(
                        [
                            wavelength[0].value,
                            (wavelength[n_wavelength - 2].value * 2.0 + wavelength[n_wavelength - 1].value) / 3.0,
                        ]
                    )  # 2 wavelengths in nm
                    # Calculate baseline SG area for scaling purposes
                    area_sg = iris_response["GEOM_AREA"]
                    for n in range(len(iris_response["INDEX_EL_SG"][nuv, :])):
                        index_values2 = iris_response["INDEX_EL_SG"][nuv, n]
                        area_sg = area_sg * iris_response["ELEMENTS"][index_values2].trans
                    # SG and SJI areas at wavelength
                    interpol_sg = scipy.interpolate.interp1d(
                        iris_response["LAMBDA"],
                        np.squeeze(area_sg),
                        fill_value="extrapolate",
                    )
                    area_sg2 = interpol_sg(wavelength)
                    area_sj2 = np.zeros((2, 2))
                    for n in range(2):
                        interpol_sji = scipy.interpolate.interp1d(
                            iris_response["LAMBDA"],
                            np.squeeze(area_sji[n]),
                            fill_value="extrapolate",
                        )
                        area_sj2[n, :] = interpol_sji(wavelength)
                    # Calculate the normalized slant function scal, apply to asji
                    for n in range(n_time_obs):
                        # Best-estimate slant, i.e., eff.area @ wavelength / baseline SG @ wavelength
                        interpol_sg2 = scipy.interpolate.interp1d(
                            iris_response["LAMBDA"],
                            np.squeeze(iris_response["AREA_SG"][0, :]),
                            fill_value="extrapolate",
                        )
                        sca2 = interpol_sg2(wavelength) / area_sg2
                        # Normalize slant so that total(wei*asj2*sca2)/total(wei*asj2)=1
                        for m in range(2):
                            sca2n = sca2 * np.sum(weight * area_sj2[m, :]) / np.sum(weight * area_sj2[m, :] * sca2)
                            interpol_sca = scipy.interpolate.interp1d(
                                wavelength, np.squeeze(sca2n), fill_value="extrapolate"
                            )
                            sca1n = interpol_sca(iris_response["LAMBDA"])
                            sca1n = np.clip(sca1n, a_min=0, a_max=None)
                            iris_response["AREA_SJI"][m] = area_sji[m] * sca1n
                else:
                    # NUV: essentially same calculation as r.version=3
                    for n in range(n_time_obs):
                        iris_response["AREA_SJI"] = [Quantity(x, unit=u.cm**2) for x in iris_response["AREA_SJI"]]
                        area_sji = [x for x in area_sji]
                        iris_response["AREA_SJI"][2:4] = area_sji[:]
            for j in range(4):
                # SJI specific time dependency
                iris_fit_sji = fit_iris_xput(
                    time_obs,
                    iris_response["C_S_TIME"][j, :, :],
                    iris_response["COEFFS_SJI"][j, :, :],
                )
                for k in range(n_time_obs):
                    iris_response["AREA_SJI"][j] = iris_response["AREA_SJI"][j] * iris_fit_sji[k]

    if not isinstance(iris_response["AREA_SG"], Quantity):
        iris_response["AREA_SG"] = Quantity(iris_response["AREA_SG"], unit=u.cm**2)
    if not isinstance(iris_response["AREA_SJI"], Quantity):
        iris_response["AREA_SJI"] = Quantity(iris_response["AREA_SJI"], unit=u.cm**2)

    return iris_response


def fit_iris_xput(time_obs, time_cal_coeffs, cal_coeffs):
    """
    To calculate the coefficients of best-fit time function for throughput,
    for which there are two modes:
    1. Perform fit: supply xput and single element ``cal_coeffs``.
    2. Apply fit: supply full ``cal_coeffs``.

    The procedure involved in this function is as follows:
    1. The time difference (in years) is computed from the time_obs and time_cal_coeffs.
    2. A least-sqaures fit is performed to determine the best fit for the time-dependent
    effective areas given the time difference.

    Parameters
    ----------
    time_obs: a `numpy.array`
        A set of observation times as `astropy.time.Time` objects contained
        in a numpy array.
    time_cal_coeffs: a numpy array of floats (with exactly two columns)
        Start and end times of intervals of constant ``cal_coeffs[i]``.
    cal_coeffs: a numpy array of floats (with at least two columns)
        Coefficients of best-fit function.

    Returns
    -------
    `numpy.array`
        Yields the fit used to compute the effective area using the input times ``time_obs``.
    """
    time_obs = np.array([parse_time(time_obs, format="utime")])

    size_time_cal_coeffs = time_cal_coeffs.shape
    size_cal_coeffs = cal_coeffs.shape

    if size_time_cal_coeffs[1] != 2 or size_cal_coeffs[1] < 2:
        # Raise ValueError as time coefficient have the wrong format.
        raise ValueError("Incorrect number of elements either in time_cal_coeffs or in cal_coeffs.")

    # Some time transformations.
    # Convert the time_cal_coeffs given in the .geny file into a ``astropy.time.Time``
    # object called t_cal_coeffs, so that the time differences will be in days...
    t_cal_coeffs_flat = time_cal_coeffs.flatten()
    t_cal_coeffs = [parse_time(x, format="utime") for x in t_cal_coeffs_flat]
    t_cal_coeffs = np.array(t_cal_coeffs).reshape(size_time_cal_coeffs)

    # Exponent for transition between exp.decay intervals.
    transition_exp = 1.5

    # For loop for carrying out the least-squares fit and computation of fit output.
    fit_out = np.zeros(len(time_obs))

    for i, t in enumerate(time_obs):
        aux_cal_coeffs = np.zeros(2 * size_time_cal_coeffs[0])

        fit_out = np.zeros(len(list(time_obs)), dtype=np.float64)

        # Looking for the closest time in the calibration time intervals.
        # Differences are given in years before passing to the next stage.
        t_diff = t - t_cal_coeffs

        t_diff = t_diff.flatten()

        # To Convert to an array, quantities need to be dimensionless, hence dividing out the unit.
        t_diff = np.array([x.to(u.year).value for x in t_diff])
        w = np.where(t_diff < 0)[0][0]

        # If the t_obs is betweeen the calibration time intervals of a
        # calibration file (w % !=0) then the aux_coeffs are given by an
        # exponential (coefficient and exponential value).
        # If the t_obs is between the end calibration time interval of
        # a calibration file (cal_file_t) and the beginning calibration
        # time interval of the next calibration file (cal_file_t+1)
        # (w % 2 == 0) then, the aux_coeffs are given by 4 values
        # corresponding to a partial exponential obtained from
        # cal_file_t and a complementary exponential obtained from the
        # cal_file_t+1
        if w % 2 != 0:  # I.e., if w is not even...
            dtt_0 = 1.0
            exp_0 = np.exp(cal_coeffs[w // 2, 2] * (t_diff[w - 1]))
            aux_cal_coeffs[w - 1 : w + 1] = np.array([dtt_0, dtt_0 * exp_0])
        else:
            dtt_1 = (t_diff[w - 1] / (t_diff[w - 1] - t_diff[w])) ** transition_exp
            dtt_0 = 1.0 - dtt_1
            exp_0 = np.exp(cal_coeffs[(w // 2) - 1, 2] * (t_diff[w - 2]))
            exp_1 = np.exp(cal_coeffs[w // 2, 2] * (t_diff[w]))
            aux_cal_coeffs[w - 2 : w + 2] = np.array([dtt_0, dtt_0 * exp_0, dtt_1, dtt_1 * exp_1])
        # print(aux_cal_coeffs)
        # print(cal_coeffs[:,:2].reshape((aux_cal_coeffs.shape[0])))
        fit_out[i] = np.matmul(aux_cal_coeffs, cal_coeffs[:, :2].reshape(aux_cal_coeffs.shape[0]))

    return fit_out
