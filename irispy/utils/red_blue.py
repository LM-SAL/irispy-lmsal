"""
Red-blue asymmetry utilities for IRIS spectrogram cubes.
"""

import warnings

import numpy as np
from scipy.interpolate import make_interp_spline

import astropy.units as u
from astropy import constants
from astropy.nddata import StdDevUncertainty
from astropy.wcs import WCS

from irispy.spectrograph import RasterCollection, SpectrogramCube
from irispy.utils._spectral import (
    _QualityFlag,
    check_scaled,
    drop_extra_coords_dependent_on_axis,
    make_map_cube,
    make_spatial_template,
    standard_deviation,
)
from irispy.utils.utils import _boolean_mask, _saturated_samples

__all__ = ["RBAQualityFlag", "calculate_red_blue_asymmetry"]

# Minimum fraction of wing bins that must have finite interpolated values.
# This is a choice of this module, not a value from the literature: a pixel with fewer
# finite bins is flagged INCOMPLETE_WINGS.
_MIN_WING_COVERAGE = 0.8
# Spatial pixels whose profiles are interpolated at a time: without return_profiles only one
# chunk of interpolated profiles is held
_CHUNK_ROWS = 20000


class RBAQualityFlag(_QualityFlag):
    """
    Quality flags for the per-pixel RBA computation.
    """

    OK = (0, "ok")
    NO_FINITE_DATA = (1, "no finite data")
    PEAK_AT_EDGE = (2, "peak at spectral edge")
    TOO_FEW_POINTS = (3, "too few finite points")
    INTERP_FAILED = (4, "interpolation failed")
    PEAK_IS_ZERO = (5, "peak is zero or non-finite")
    INCOMPLETE_WINGS = (6, "incomplete red or blue wing coverage")
    LOW_SIGNAL = (7, "below min_intensity")
    SATURATED = (8, "saturated")


def _make_velocity_wcs(array_shape, velocity_axis, velocity_grid):
    """
    Return a profile WCS whose only meaningful world axis is Doppler velocity.

    Do not copy the input cube WCS here: profile cubes may have fewer axes than
    the source cube, so stale FITS-WCS axes can warn or fail during construction.
    """
    naxis = len(array_shape)
    wcs_axis = naxis - 1 - velocity_axis
    cdelt = float(np.nanmean(np.diff(velocity_grid))) if velocity_grid.size > 1 else 1.0

    wcs = WCS(naxis=naxis)
    wcs.wcs.ctype[wcs_axis] = "VELO"
    wcs.wcs.cunit[wcs_axis] = "km/s"
    wcs.wcs.crpix[wcs_axis] = 1.0
    wcs.wcs.crval[wcs_axis] = float(velocity_grid[0])
    wcs.wcs.cdelt[wcs_axis] = cdelt
    wcs.array_shape = tuple(array_shape)
    wcs.wcs.set()
    return wcs


def _make_profile_cube(cube, *, data, velocity_grid, wavelength_axis, meta, uncertainty=None, mask=None):
    return SpectrogramCube(
        data,
        wcs=_make_velocity_wcs(data.shape, wavelength_axis, velocity_grid),
        uncertainty=uncertainty,
        unit=cube.unit,
        meta=meta,
        mask=mask,
        extra_coords=drop_extra_coords_dependent_on_axis(cube.extra_coords, wavelength_axis, reindex=False),
    )


def _prepare_data(cube, wavelength_axis):
    """
    Data and errors, negatives masked, spectral axis last.
    """
    data = np.asarray(cube.data, dtype=float)
    # Negative samples are left out
    mask = _boolean_mask(cube)
    data = np.where((data < 0) | mask, np.nan, data)

    errors = standard_deviation(cube)
    if errors is not None:
        errors = np.where(mask, np.nan, errors)

    data = np.moveaxis(data, wavelength_axis, -1)
    if errors is not None:
        errors = np.moveaxis(errors, wavelength_axis, -1)
    return data, errors


def calculate_red_blue_asymmetry(
    cube,
    *,
    rest_wavelength=None,
    velocity_range=(50, 150) * u.km / u.s,
    dv=10 * u.km / u.s,
    spline_degree=3,
    min_intensity=None,
    return_profiles=True,
):
    """
    Calculate red-blue asymmetry maps from a spectrogram cube.

    The asymmetry is computed for each spatial pixel as
    :math:`(I_R - I_B) / I_p`, where ``I_R`` and ``I_B`` are mean intensities
    in matching red and blue velocity ranges and ``I_p`` is the interpolated
    peak intensity.

    Parameters
    ----------
    cube : `irispy.spectrograph.SpectrogramCube`
        Input spectrogram cube.
    rest_wavelength : `astropy.units.Quantity`, optional
        Rest wavelength used to convert the spectral axis to Doppler velocity.
        If omitted, read from ``cube.meta.rest_wavelength`` (requires an
        `~irispy.meta.SGMeta` instance with a ``TWAVE<n>`` FITS keyword).
    velocity_range : `astropy.units.Quantity`, optional
        Two positive velocities. The red wing spans them and the blue wing their negatives,
        both measured from the brightest sample of each profile.
    dv : `astropy.units.Quantity`, optional
        Velocity spacing for the interpolated profile.
    spline_degree : `int`, optional
        Spline degree for `scipy.interpolate.make_interp_spline`.
    min_intensity : `float` or `astropy.units.Quantity`, optional
        Minimum peak intensity required for a pixel to be processed.
        Pixels below this threshold are assigned quality flag
        `~irispy.utils.red_blue.RBAQualityFlag.LOW_SIGNAL`.
    return_profiles : `bool`, optional
        If `True`, include 3D ``"observed_profile"`` and
        ``"interpolated_profile"`` cubes in the output.

    Returns
    -------
    `irispy.spectrograph.RasterCollection`
        Always contains ``"red_blue_asymmetry"`` and ``"quality"`` 2D maps.
        ``"red_blue_asymmetry_error"`` is added when the input cube has uncertainty.
        ``"observed_profile"`` and ``"interpolated_profile"`` 3D cubes are added
        when ``return_profiles=True``. The interpolated profile velocity axis is
        relative to the brightest sample of each profile; the observed profile velocity
        axis is relative to ``rest_wavelength``.

    Raises
    ------
    ValueError
        If no velocity of the interpolated profile, which is spaced by ``dv``, falls in
        ``velocity_range``, as when the range is narrower than ``dv``.

    Notes
    -----
    Red-blue asymmetries were used to find the upflows of the transition region and corona by
    :cite:t:`depontieu2009` and :cite:t:`tian2011`. This function is not a port of either method.
    These choices are its own, and the papers should be consulted before comparing values with
    theirs: the profile is interpolated with a spline on a regular velocity grid, centered on its
    brightest sample and not on a fitted centroid or the rest wavelength, the
    asymmetry is the difference of the mean intensities of the two wings over the interpolated
    peak intensity, and the default wing range of 50 to 150 km/s and the requirement that
    80% of the wing bins are finite have no source.

    Pixels with a saturated sample, as the mask and ``cube.meta["saturated"]`` mark the samples the
    readers find at the level 2 ceiling, are flagged `~irispy.utils.red_blue.RBAQualityFlag.SATURATED`.
    Subtract a background first with `~irispy.utils.spectrograph.subtract_background` if the line
    sits on one.
    """
    check_scaled(cube)
    if rest_wavelength is None:
        rest_wavelength = getattr(cube.meta, "rest_wavelength", None)
        if rest_wavelength is None:
            msg = (
                "rest_wavelength was not provided and could not be read from cube.meta. "
                "Pass rest_wavelength explicitly or provide a cube with SGMeta containing TWAVE keywords."
            )
            raise ValueError(msg)
    rest_wavelength = u.Quantity(rest_wavelength).to(u.nm)
    if rest_wavelength.shape != () or not np.isfinite(rest_wavelength.value) or rest_wavelength <= 0 * u.nm:
        msg = "rest_wavelength must be a positive finite scalar wavelength"
        raise ValueError(msg)

    velocity_range = u.Quantity(velocity_range).to(u.km / u.s)
    if velocity_range.shape != (2,):
        msg = "velocity_range must contain two velocities"
        raise ValueError(msg)
    velocity_low, velocity_high = velocity_range
    if velocity_low < 0 or velocity_high <= velocity_low:
        msg = "velocity_range must be positive and increasing"
        raise ValueError(msg)

    dv = u.Quantity(dv).to_value(u.km / u.s)
    if dv <= 0:
        msg = "dv must be positive"
        raise ValueError(msg)
    spline_degree = int(spline_degree)
    if spline_degree < 0:
        msg = "spline_degree must be a non-negative integer"
        raise ValueError(msg)

    wavelength_axis = cube.wavelength_axis
    wavelengths = cube.axis_world_coords(wavelength_axis)[0].to(u.nm)
    rest_wavelength = rest_wavelength.to(wavelengths.unit)
    velocity = ((wavelengths - rest_wavelength) / rest_wavelength * constants.c).to_value(u.km / u.s)
    interp_extent = velocity_high.to_value(u.km / u.s) + dv
    interp_velocity = np.arange(-interp_extent, interp_extent + dv, dv)
    velocity_range_kms = (velocity_low.to_value(u.km / u.s), velocity_high.to_value(u.km / u.s))
    low, high = velocity_range_kms
    red_mask = (interp_velocity >= low) & (interp_velocity <= high)
    blue_mask = (interp_velocity >= -high) & (interp_velocity <= -low)
    if not (red_mask.any() and blue_mask.any()):
        msg = f"No velocity of the interpolated profile, spaced by dv={dv} km/s, falls in velocity_range=({low}, {high}) km/s"
        raise ValueError(msg)

    data, errors = _prepare_data(cube, wavelength_axis)
    output_shape = data.shape[:-1]
    # One row per spatial pixel
    profiles = data.reshape(-1, data.shape[-1])
    profile_errors = None if errors is None else errors.reshape(profiles.shape)

    quality = np.full(len(profiles), RBAQualityFlag.OK, dtype=np.uint8)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        raw_peak = np.nanmax(profiles, axis=-1)
    if min_intensity is not None:
        quality[raw_peak < u.Quantity(min_intensity, cube.unit).value] = RBAQualityFlag.LOW_SIGNAL
    quality[_saturated_samples(cube).any(axis=wavelength_axis).ravel()] = RBAQualityFlag.SATURATED

    min_points = spline_degree + 1
    finite = np.isfinite(profiles)
    finite_velocity = np.isfinite(velocity)
    todo = quality == RBAQualityFlag.OK
    quality[todo & ~finite.any(axis=-1)] = RBAQualityFlag.NO_FINITE_DATA
    todo &= quality == RBAQualityFlag.OK
    peak_index = np.argmax(np.where(finite, profiles, -np.inf), axis=-1)
    quality[todo & ((peak_index == 0) | (peak_index == profiles.shape[1] - 1))] = RBAQualityFlag.PEAK_AT_EDGE
    todo &= quality == RBAQualityFlag.OK
    # The spline goes through the finite samples at finite velocities from the peak
    points = finite & finite_velocity
    n_points = np.where(finite_velocity[peak_index], points.sum(axis=-1), 0)
    quality[todo & (n_points < min_points)] = RBAQualityFlag.TOO_FEW_POINTS
    todo &= quality == RBAQualityFlag.OK

    def interpolate(x, y):
        return make_interp_spline(x, y, k=spline_degree, axis=-1)(interp_velocity, extrapolate=False)

    def asymmetries(profiles, profile_errors, peak_index, points, todo, quality):
        """
        The asymmetries and errors of a chunk of rows, flagging ``quality`` in place,
        with their interpolated profiles and errors.
        """
        interpolated = np.full((len(profiles), interp_velocity.size), np.nan)
        interpolated_errors = None if profile_errors is None else np.full_like(interpolated, np.nan)
        # Profiles without gaps share their velocity grid, so one spline fits all of those whose peak
        # is at the same sample; a gap leaves a profile its own grid, and its own spline.
        complete = todo & points.all(axis=-1)
        if profile_errors is not None:
            complete &= np.isfinite(profile_errors).all(axis=-1)
        order = np.argsort(velocity)
        for index in np.unique(peak_index[complete]):
            rows = np.flatnonzero(complete & (peak_index == index))
            x = velocity[order] - velocity[index]
            try:
                interpolated[rows] = interpolate(x, profiles[rows][:, order])
                if profile_errors is not None:
                    interpolated_errors[rows] = interpolate(x, profile_errors[rows][:, order])
            except ValueError:
                quality[rows] = RBAQualityFlag.INTERP_FAILED
                interpolated[rows] = np.nan
        for row in np.flatnonzero(todo & ~complete):
            kept = points[row]
            shifted_velocity = velocity[kept] - velocity[peak_index[row]]
            order = np.argsort(shifted_velocity)
            x = shifted_velocity[order]
            samples = profiles[row][kept][order]
            row_errors = None if profile_errors is None else profile_errors[row][kept][order]
            finite_errors = None if row_errors is None else np.isfinite(row_errors)
            try:
                if row_errors is not None and finite_errors.all():
                    # Errors on the same grid: one spline through both
                    interpolated[row], interpolated_errors[row] = interpolate(x, [samples, row_errors])
                else:
                    interpolated[row] = interpolate(x, samples)
                    if row_errors is not None and finite_errors.sum() >= min_points:
                        interpolated_errors[row] = interpolate(x[finite_errors], row_errors[finite_errors])
            except ValueError:
                quality[row] = RBAQualityFlag.INTERP_FAILED
                interpolated[row] = np.nan
        if interpolated_errors is not None:
            interpolated_errors = np.where(interpolated_errors >= 0, interpolated_errors, np.nan)

        finite_interpolated = np.isfinite(interpolated)
        incomplete = ((finite_interpolated & red_mask).sum(axis=-1) < _MIN_WING_COVERAGE * red_mask.sum()) | (
            (finite_interpolated & blue_mask).sum(axis=-1) < _MIN_WING_COVERAGE * blue_mask.sum()
        )
        quality[(quality == RBAQualityFlag.OK) & incomplete] = RBAQualityFlag.INCOMPLETE_WINGS
        rows = np.flatnonzero(quality == RBAQualityFlag.OK)  # with finite samples in both wings
        profile = interpolated[rows]
        red_intensity = np.nanmean(profile[:, red_mask], axis=-1)
        blue_intensity = np.nanmean(profile[:, blue_mask], axis=-1)
        peak = np.nanmax(profile, axis=-1)
        zero = ~np.isfinite(peak) | np.isclose(peak, 0)
        quality[rows[zero]] = RBAQualityFlag.PEAK_IS_ZERO
        rows, profile, peak = rows[~zero], profile[~zero], peak[~zero]
        numerator = (red_intensity - blue_intensity)[~zero]
        red_blue = np.full(len(profiles), np.nan)
        red_blue[rows] = numerator / peak
        red_blue_error = np.full(len(profiles), np.nan)
        if interpolated_errors is not None:
            error = interpolated_errors[rows]

            def wing_error(wing):
                count = np.isfinite(error[:, wing]).sum(axis=-1)
                total = np.sqrt(np.nansum(error[:, wing] ** 2, axis=-1))
                return np.divide(total, count, out=np.full(count.shape, np.nan), where=count > 0)

            numerator_error = np.sqrt(wing_error(red_mask) ** 2 + wing_error(blue_mask) ** 2)
            peak_error = error[np.arange(len(rows)), np.nanargmax(profile, axis=-1)]
            numerator_is_zero = np.isclose(numerator, 0)
            valid = np.isfinite(numerator_error) & (numerator_is_zero | np.isfinite(peak_error))
            variance = (numerator_error / peak) ** 2
            variance = variance + np.where(numerator_is_zero, 0, (numerator * peak_error / peak**2) ** 2)
            red_blue_error[rows[valid]] = np.sqrt(variance[valid])
        return red_blue, red_blue_error, interpolated, interpolated_errors

    red_blue = np.full(len(profiles), np.nan)
    red_blue_error = np.full(len(profiles), np.nan)
    interpolated = interpolated_errors = None
    if return_profiles:
        interpolated = np.full((len(profiles), interp_velocity.size), np.nan)
        interpolated_errors = None if profile_errors is None else np.full_like(interpolated, np.nan)
    for start in range(0, len(profiles), _CHUNK_ROWS):
        chunk = slice(start, start + _CHUNK_ROWS)
        red_blue[chunk], red_blue_error[chunk], chunk_profiles, chunk_errors = asymmetries(
            profiles[chunk],
            None if profile_errors is None else profile_errors[chunk],
            peak_index[chunk],
            points[chunk],
            todo[chunk],
            quality[chunk],
        )
        if return_profiles:
            interpolated[chunk] = chunk_profiles
            if interpolated_errors is not None:
                interpolated_errors[chunk] = chunk_errors

    red_blue = red_blue.reshape(output_shape)
    red_blue_error = red_blue_error.reshape(output_shape)
    quality = quality.reshape(output_shape)

    meta = {
        "rba_rest_wavelength": rest_wavelength.to_value(u.nm),
        "rba_rest_wavelength_unit": "nm",
        "rba_velocity_range": velocity_range_kms,
        "rba_dv": dv,
        "rba_spline_degree": spline_degree,
    }

    template = make_spatial_template(cube, wavelength_axis)

    def _make_cube(values, unit, *, mask=None):
        c = make_map_cube(template, values, unit, mask=mask)
        c.meta.update(meta)
        return c

    cubes = [
        ("red_blue_asymmetry", red_blue, u.dimensionless_unscaled),
        ("quality", quality, u.dimensionless_unscaled),
    ]
    if errors is not None:
        cubes.append(("red_blue_asymmetry_error", red_blue_error, u.dimensionless_unscaled))

    result_cubes = []
    for name, values, unit in cubes:
        mask = None if name == "quality" else ~np.isfinite(values)
        result_cubes.append((name, _make_cube(values, unit, mask=mask)))

    if return_profiles:
        observed_data = np.moveaxis(data, -1, wavelength_axis)
        interpolated_data = np.moveaxis(interpolated.reshape(*output_shape, -1), -1, wavelength_axis)
        observed_mask = ~np.isfinite(observed_data)
        interpolated_mask = ~np.isfinite(interpolated_data)
        observed_uncertainty = (
            StdDevUncertainty(np.moveaxis(errors, -1, wavelength_axis)) if errors is not None else None
        )
        interpolated_uncertainty = (
            StdDevUncertainty(np.moveaxis(interpolated_errors.reshape(*output_shape, -1), -1, wavelength_axis))
            if interpolated_errors is not None
            else None
        )

        result_cubes.extend(
            [
                (
                    "observed_profile",
                    _make_profile_cube(
                        cube,
                        data=observed_data,
                        velocity_grid=velocity,
                        wavelength_axis=wavelength_axis,
                        meta={**cube.meta, **meta, "rba_profile": "observed"},
                        uncertainty=observed_uncertainty,
                        mask=observed_mask,
                    ),
                ),
                (
                    "interpolated_profile",
                    _make_profile_cube(
                        cube,
                        data=interpolated_data,
                        velocity_grid=interp_velocity,
                        wavelength_axis=wavelength_axis,
                        meta={**cube.meta, **meta, "rba_profile": "interpolated_peak_centered"},
                        uncertainty=interpolated_uncertainty,
                        mask=interpolated_mask,
                    ),
                ),
            ]
        )

    return RasterCollection(result_cubes)
