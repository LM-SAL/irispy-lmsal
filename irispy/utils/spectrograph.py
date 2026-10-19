"""
This module provides general utility functions called by code in spectrograph.
"""

import numpy as np

import astropy.units as u
from astropy import constants
from astropy.nddata import StdDevUncertainty

from irispy.sji import SJICube
from irispy.spectrograph import SpectrogramCube, SpectrogramCubeSequence
from irispy.utils._spectral import _fit_background, check_scaled, make_map_cube, standard_deviation
from irispy.utils.constants import RADIANCE_UNIT, RADIANCE_UNIT_PER_HZ
from irispy.utils.response import get_interpolated_effective_area, get_latest_response
from irispy.utils.utils import _saturated_samples

__all__ = [
    "calculate_dn_to_radiance_factor",
    "convert_photons_per_sec_to_radiance",
    "radiation_temperature",
    "radiometric_calibration",
    "subtract_background",
]


def radiometric_calibration(
    cube: SpectrogramCube | SpectrogramCubeSequence | SJICube,
    *,
    line=None,
) -> SpectrogramCube | SpectrogramCubeSequence | SJICube:
    """
    Performs radiometric calibration on the input cube or cube sequence.

    This takes into consideration also the observation time and uses the latest response.

    The data is also exposure time corrected during the conversion.

    This takes into account the spectral dispersion and solid angle of the pixels
    based on the WCS, which is different from the IDL code that does not take spectral
    dispersion into account. If you want the same results as the IDL code, you can multiply
    the output by the spectral dispersion.

    Parameters
    ----------
    cube : `irispy.spectrograph.SpectrogramCube` | `irispy.spectrograph.SpectrogramCubeSequence` | `irispy.sji.SJICube`
        The input cube to be calibrated.
    line : `astropy.units.Quantity`, optional
        Slit-jaw cubes only: the wavelength of the one line all the counts come from, which gives the
        line intensity in :math:`erg s^{-1} sr^{-1} cm^{-2}`. Defaults to `None`, a continuum flat across
        the passband, which gives :math:`erg s^{-1} sr^{-1} cm^{-2} Å^{-1}`.
        :ref:`irispy-tutorial-calibration` says when each holds.

    Returns
    -------
    `irispy.spectrograph.SpectrogramCube`, `irispy.spectrograph.SpectrogramCubeSequence` or `irispy.sji.SJICube`
        New cube in new units.

    Raises
    ------
    ValueError
        If ``line`` is given for a spectrograph cube or lies where the slit-jaw passband has no
        effective area, or if the cube is not from an IRIS slit-jaw passband, such as an AIA cube.

    Notes
    -----
    This is designed to do the same as `iris2/iris_calib_spectrum.pro <https://hesperia.gsfc.nasa.gov/ssw/iris/idl/lmsal/iris2/iris_calib_spectrum.pro>`__ IDL code.

    The calibration output has been confirmed to provide the same results as those provided
    by the SolarSoft IDL routine `IRIS_CALIB <https://hesperia.gsfc.nasa.gov/ssw/iris/idl/nrl/iris_calib.pro>`__.
    The major difference being that the output here is accounting for the wavelength, which is why the units
    here are :math:`erg s^{-1} sr^{-1} cm^{-2} Å^{-1}` and not :math:`erg s^{-1} sr^{-1} cm^{-2}`.
    Notice the extra :math:`Å^{-1}` in the units.

    The response file only defines the effective area within the nominal spectral
    ranges of each band, while some observations read out wider spectral windows
    (e.g., a full-CCD Si IV window reaching blueward of ~1389 Å). The upstream IDL
    code divides by the zero effective area there and silently returns infinities;
    here those wavelengths are instead set to NaN in the output data. The mask is
    left untouched (it is copied from the input cube) so nothing is hidden beyond
    what was already flagged as bad in the level 2 data.
    """
    if isinstance(cube, SpectrogramCubeSequence):
        return SpectrogramCubeSequence([radiometric_calibration(c, line=line) for c in cube])
    check_scaled(cube)
    if isinstance(cube, SJICube):
        return _calibrate_sji(cube, line)
    if line is not None:
        msg = "line is for slit-jaw cubes; spectra are calibrated at each of their wavelengths"
        raise ValueError(msg)
    detector_type = cube.meta.detector_band
    spectral_dispersion_per_pixel = cube.spectral_dispersion
    solid_angle = cube.solid_angle
    # Get wavelength for each pixel.
    wavelength_axis_index = cube.wavelength_axis
    wavelength = cube.axis_world_coords(wavelength_axis_index)[0]
    time_obs = cube.meta.date_reference
    iris_response = get_latest_response(time_obs)
    exp_corrected_cube = cube.apply_exposure_time_correction()
    # Convert to radiance units.
    data_quantities = (exp_corrected_cube.data * exp_corrected_cube.unit.to(u.photon / u.s) * (u.photon / u.s),)
    if exp_corrected_cube.uncertainty is not None:
        uncertainty = (
            exp_corrected_cube.uncertainty.array * exp_corrected_cube.unit.to(u.photon / u.s) * (u.photon / u.s)
        )
        data_quantities += (uncertainty,)
    new_data_quantities = convert_photons_per_sec_to_radiance(
        data_quantities=data_quantities,
        iris_response=iris_response,
        wavelength=wavelength,
        detector_type=detector_type,
        spectral_dispersion_per_pixel=spectral_dispersion_per_pixel,
        solid_angle=solid_angle,
    )
    new_data = new_data_quantities[0].value
    new_uncertainty = StdDevUncertainty(new_data_quantities[1].value) if len(new_data_quantities) > 1 else None
    return make_map_cube(cube, new_data, new_data_quantities[0].unit, uncertainty=new_uncertainty)


def _calibrate_sji(cube, line):
    response = get_latest_response(cube.meta.date_reference)
    # NAME_SJI holds the first three digits of each passband, "133" to "283"
    names = [name.decode() for name in response["NAME_SJI"]]
    if cube.meta.spectral_window[:3] not in names:
        msg = f"There is no IRIS slit-jaw response for {cube.meta.spectral_window}"
        raise ValueError(msg)
    band = names.index(cube.meta.spectral_window[:3])
    area, wavelength = response["AREA_SJI"][band], response["LAMBDA"]
    # The intensity times t Ω, per photon counted
    if line is None:
        # DN2PHOT DN = t Ω ∫ A(λ) I_λ λ / (hc) dλ, solved for a constant I_λ
        per_photon = constants.h * constants.c / np.trapezoid(area * wavelength, wavelength)
        unit = RADIANCE_UNIT
    else:
        # The response table is zero or absent where the passband does not reach
        line_area = np.interp(line, wavelength, area, left=0, right=0)
        if not line_area > 0:
            msg = f"The {cube.meta.spectral_window} slit-jaw passband has no effective area at line={line}"
            raise ValueError(msg)
        per_photon = constants.h * constants.c / line / line_area
        unit = RADIANCE_UNIT * u.AA
    solid_angle = cube.meta["CDELT1"] * cube.meta["CDELT2"] * u.arcsec**2
    # The cube's unit, DN_IRIS_SJI_<passband>, holds the passband's photons per DN, those of DN2PHOT_SJI
    return (cube.apply_exposure_time_correction() * (per_photon / u.photon / solid_angle)).to(unit)


def radiation_temperature(
    cube: SpectrogramCube | SpectrogramCubeSequence,
) -> SpectrogramCube | SpectrogramCubeSequence:
    r"""
    Converts a radiometrically calibrated cube or cube sequence to radiation temperature.

    The radiation (or brightness) temperature is the temperature of the blackbody whose
    Planck function equals the observed specific intensity, :math:`I_\nu = B_\nu(T_\mathrm{rad})`
    (:cite:t:`rybicki1985`, Section 1.5). Inverting the Planck function at each wavelength gives

    .. math::

       T_\mathrm{rad} = \frac{h \nu / k}{\ln \left( 1 + 2 h \nu^3 / (c^2 I_\nu) \right)}.

    :cite:t:`leenaarts2013` and :cite:t:`pereira2013` express Mg II h & k and NUV
    intensities this way.

    Parameters
    ----------
    cube : `irispy.spectrograph.SpectrogramCube` | `irispy.spectrograph.SpectrogramCubeSequence`
        Cube in radiance per unit wavelength or frequency, e.g. the output of
        `~irispy.utils.spectrograph.radiometric_calibration`.

    Returns
    -------
    `irispy.spectrograph.SpectrogramCube` or `irispy.spectrograph.SpectrogramCubeSequence`
        New cube in K, with a `~astropy.nddata.StdDevUncertainty` if ``cube`` has an uncertainty.

    Notes
    -----
    * Do not use the `astropy.units.brightness_temperature` equivalency for this: it is the
      Rayleigh-Jeans limit, which is wrong by orders of magnitude in the UV.
    * The radiation temperature is not the gas temperature. They are equal only for optically thick
      radiation whose source function is the Planck function (:cite:t:`rybicki1985`, Sections 1.4
      and 1.5), which holds only approximately for Mg II h & k :cite:p:`leenaarts2013` and not at
      all for optically thin lines.
    * Samples with zero, negative or NaN radiance are NaN and masked, and +Inf samples are +Inf and
      masked; the input mask is kept.
      Noise just above zero still gives several thousand kelvin.
    * The uncertainty is propagated to first order, so it is unreliable where it is comparable to the radiance.
    """
    if isinstance(cube, SpectrogramCubeSequence):
        return SpectrogramCubeSequence([radiation_temperature(c) for c in cube])
    check_scaled(cube)
    wavelength = reshape_1d_wavelength_dimensions_for_broadcast(
        cube.axis_world_coords(cube.wavelength_axis)[0], cube.data.ndim
    )
    try:
        to_radiance = cube.unit.to(RADIANCE_UNIT_PER_HZ, equivalencies=u.spectral_density(wavelength))
    except u.UnitConversionError:
        msg = (
            f"The cube must be in radiance per unit wavelength or frequency, not {cube.unit}; "
            "convert it with radiometric_calibration first"
        )
        raise ValueError(msg) from None
    frequency = wavelength.to(u.Hz, equivalencies=u.spectral())
    temperature_scale = (constants.h * frequency / constants.k_B).to_value(u.K)
    planck_scale = (2 * constants.h * frequency**3 / constants.c**2 / u.sr).to_value(RADIANCE_UNIT_PER_HZ)
    radiance = np.asarray(cube.data, dtype=float) * to_radiance
    sigma = standard_deviation(cube)
    uncertainty = None
    with np.errstate(divide="ignore", invalid="ignore"):
        temperature = np.where(radiance > 0, temperature_scale / np.log1p(planck_scale / radiance), np.nan)
        if sigma is not None:
            # First-order propagation with dT/dI from differentiating the expression above.
            derivative = temperature**2 * planck_scale / (temperature_scale * radiance * (radiance + planck_scale))
            uncertainty = StdDevUncertainty(derivative * sigma * to_radiance)
    return make_map_cube(cube, temperature, u.K, mask_invalid=True, uncertainty=uncertainty)


def convert_photons_per_sec_to_radiance(
    *,
    data_quantities,
    iris_response,
    wavelength,
    detector_type,
    spectral_dispersion_per_pixel,
    solid_angle,
):
    """
    Converts data quantities from counts/s to radiance.

    Parameters
    ----------
    data_quantities: iterable of `astropy.units.Quantity`
        Quantities to be converted.  Must have units of counts/s or
        radiance equivalent counts, e.g. erg / cm**2 / s / sr / Angstrom.
    iris_response: dict
        The IRIS response data loaded from `irispy.utils.response.get_latest_response`.
    wavelength: `astropy.units.Quantity`
        Wavelength at each element along spectral axis of data quantities.
    detector_type: `str`
        Detector type: 'FUV', 'NUV', or 'SJI'.
    spectral_dispersion_per_pixel: scalar `astropy.units.Quantity`
        Spectral dispersion (wavelength width) of a pixel.
    solid_angle: scalar `astropy.units.Quantity`
        Solid angle corresponding to a pixel.

    Returns
    -------
    `list` of `astropy.units.Quantity`
        Data quantities converted to radiance.

    Notes
    -----
    This is designed to do the same as `nrl/iris_calib.pro <https://hesperia.gsfc.nasa.gov/ssw/iris/idl/nrl/iris_calib.pro>`__ IDL code.
    The difference is that this function takes into account the spectral dispersion which the IDL code
    does not.
    To get the same results as the IDL code, you can multiply the output by the spectral dispersion
    or set the keyword to have the value of 1 Angstrom.
    """
    for i, data in enumerate(data_quantities):
        if data.unit != u.photon / u.s:
            msg = (
                f"Invalid unit provided. Unit must be equivalent to {u.photon / u.s}. "
                f"Error found for {i}th element of ``data_quantities`` with unit: {data.unit}"
            )
            raise ValueError(
                msg,
            )
    photons_per_sec_to_radiance_factor = calculate_dn_to_radiance_factor(
        iris_response=iris_response,
        wavelength=wavelength,
        detector_type=detector_type,
        spectral_dispersion_per_pixel=spectral_dispersion_per_pixel,
        solid_angle=solid_angle,
    )
    # Change shape of arrays so they are compatible for broadcasting
    # with data and uncertainty arrays.
    photons_per_sec_to_radiance_factor = reshape_1d_wavelength_dimensions_for_broadcast(
        photons_per_sec_to_radiance_factor,
        data_quantities[0].ndim,
    )
    return [(data * photons_per_sec_to_radiance_factor).to(RADIANCE_UNIT) for data in data_quantities]


def calculate_dn_to_radiance_factor(
    *,
    iris_response,
    wavelength,
    detector_type,
    spectral_dispersion_per_pixel,
    solid_angle,
):
    """
    Calculates multiplicative factor that converts counts/s to radiance for given
    wavelengths.

    Parameters
    ----------
    iris_response: dict
        The IRIS response data loaded from `irispy.utils.response.get_latest_response`.
    wavelength: `astropy.units.Quantity`
        Wavelengths for which counts/s-to-radiance factor is to be calculated
    detector_type: `str`
        Detector type: 'FUV' or 'NUV'.
    spectral_dispersion_per_pixel: scalar `astropy.units.Quantity`
        Spectral dispersion (wavelength width) of a pixel.
    solid_angle: scalar `astropy.units.Quantity`
        Solid angle corresponding to a pixel.

    Returns
    -------
    `astropy.units.Quantity`
        Multiplicative conversion factor from counts/s to radiance units
        for input wavelengths.

    Notes
    -----
    The term "multiplicative" refers to the fact that the conversion factor calculated by the
    `.calculate_dn_to_radiance_factor` function is used to multiply the counts per
    second (cps) data to obtain the radiance data. In other words, the conversion factor is a
    scaling factor that is applied to the cps data to convert it to radiance units.
    """
    # Get effective area and interpolate to observed wavelength grid.
    eff_area_interp = get_interpolated_effective_area(
        iris_response,
        detector_type,
        wavelength,
    )
    # Return radiometric converted data assuming input data is in units of photons/s.
    return (
        constants.h
        * constants.c
        / wavelength
        / u.photon
        / spectral_dispersion_per_pixel
        / eff_area_interp
        / solid_angle
    )


def reshape_1d_wavelength_dimensions_for_broadcast(wavelength, n_data_dim):
    if n_data_dim == 1:
        pass
    elif n_data_dim == 2:
        wavelength = wavelength[np.newaxis, :]
    elif n_data_dim == 3:
        wavelength = wavelength[np.newaxis, np.newaxis, :]
    else:
        msg = "IRISSpectrogram dimensions must be 2 or 3."
        raise ValueError(msg)
    return wavelength


def subtract_background(cube, windows, *, degree=1):
    r"""
    Fit a polynomial to line-free wavelength windows of each spectrum and subtract it.

    Only unmasked, finite samples within ``windows`` determine the polynomial.
    The fitted background is then evaluated and subtracted at every wavelength in the spectrum.

    Parameters
    ----------
    cube : `irispy.spectrograph.SpectrogramCube`
        Input cube with a wavelength axis.
    windows : `astropy.units.Quantity`
        A ``(lower, upper)`` wavelength window, or an ``(n, 2)`` array of them, holding no line.
    degree : `int`, optional
        Degree of the polynomial in wavelength, fitted to the unmasked, finite samples in the
        windows by unweighted least squares. Defaults to 1, a straight line.

    Returns
    -------
    `irispy.spectrograph.SpectrogramCube`
        ``cube`` minus the fitted background, with the same mask, uncertainty and coordinates.
        Spectra with no more than ``degree`` samples to fit are NaN, except their saturated samples.

    Notes
    -----
    The error of the fitted background is shared by all the samples of a spectrum, so it adds
    up coherently in a sum over wavelength, which a
    `~astropy.nddata.StdDevUncertainty` of independent samples cannot hold. For a constant
    (``degree=0``) fitted to :math:`n` samples of error :math:`\sigma`, it is :math:`\sigma / \sqrt{n}`.

    Saturated samples (see `~irispy.utils.moments.calculate_moments`) are masked, so left out of the
    fit, and keep their value in spectra with too few samples to fit, so that
    `~irispy.utils.moments.calculate_moments` still flags a spectrum saturated across its line and its
    windows.
    """
    check_scaled(cube)
    background = _fit_background(cube, windows, degree)
    # Their value minus NaN would be NaN
    background[np.isnan(background) & _saturated_samples(cube)] = 0
    return cube - u.Quantity(background, cube.unit)
