"""
Helpers for IRIS density diagnostics.
"""

import numpy as np

import astropy.units as u
from astropy.nddata import NDDataRef, StdDevUncertainty

from irispy.utils._spectral import standard_deviation

__all__ = ["density_diagnostic"]


def _map_ratio_to_quantity(observed_ratio, quantity, theoretical_ratio, fill_value=np.nan):
    """
    ``quantity`` at ``observed_ratio`` on the monotonic ``theoretical_ratio``,
    interpolated linearly, and ``fill_value``, one value or a ``(left, right)`` pair,
    off it.
    """
    if np.all(np.diff(theoretical_ratio) <= 0.0):
        quantity = quantity[::-1]
        theoretical_ratio = theoretical_ratio[::-1]
    order = np.argsort(theoretical_ratio)
    theoretical_ratio, unique_index = np.unique(theoretical_ratio[order], return_index=True)
    quantity = quantity[order][unique_index]
    left, right = fill_value if isinstance(fill_value, tuple) else (fill_value, fill_value)
    return u.Quantity(
        np.interp(observed_ratio, theoretical_ratio, quantity.value, left=left, right=right), quantity.unit
    )


def _with_standard_deviation(intensity):
    # An UnknownUncertainty, which astropy cannot propagate, is taken to be a standard deviation
    intensity = NDDataRef(intensity)
    if (sigma := standard_deviation(intensity)) is None:
        return intensity
    return NDDataRef(intensity.data, uncertainty=StdDevUncertainty(sigma), mask=intensity.mask, unit=intensity.unit)


def density_diagnostic(intensity_numerator, intensity_denominator, density_grid, theoretical_ratio):
    """
    Map an observed IRIS line ratio onto a theoretical density curve.

    Parameters
    ----------
    intensity_numerator, intensity_denominator : array-like, `~astropy.units.Quantity` or `~astropy.nddata.NDData`
        Measured integrated intensities of the numerator and denominator lines on the same pixels, such as two
        `~ndcube.NDCube` maps of one raster. The uncertainties of `~astropy.nddata.NDData` intensities give the
        density bounds, and their masked pixels give NaN. An intensity without uncertainty counts as exact.
    density_grid : `~astropy.units.Quantity`
        Densities at which ``theoretical_ratio`` is sampled. It can be any quantity the ratio depends
        on, such as the temperature, as in SolarSoft's
        `iris_te.pro <https://sohoftp.nascom.nasa.gov/solarsoft/iris/idl/sao/util/tian/iris_te.pro>`__;
        the ``"density"`` results are then values of that quantity.
    theoretical_ratio : array-like or `~astropy.units.Quantity`
        Theoretical ratio of the numerator line to the denominator line at each ``density_grid``, such as one
        computed with `fiasco.Ion.contribution_function`, as in the
        :ref:`O IV density example <sphx_glr_generated_gallery_analysis_06_oiv_density_diagnostics.py>`, or the
        O IV 1399.78 / 1401.16 Å table in SolarSoft's ``iris_ne_oiv.pro``.
        Where the curve turns, the branch that holds the most observed ratios is used.

    Returns
    -------
    `dict`
        * ``"ratio"`` — the observed ratio, NaN where the denominator is 0 or a pixel is masked
        * ``"ratio_uncertainty"`` — its uncertainty, NaN where the ratio is NaN, `None` for intensities without
          uncertainties
        * ``"density"`` — the density at the ratio, NaN off the curve
        * ``"density_lower"``, ``"density_upper"`` — the densities at the ratio minus and plus its uncertainty,
          `None` for intensities without uncertainties. A bound at an end of the curve means the error bar runs
          off the curve.
        * ``"density_grid"``, ``"theoretical_ratio"`` — the curve used: the finite samples of its branch, which
          can be shorter than ``density_grid``

    Raises
    ------
    ValueError
        If ``density_grid`` and ``theoretical_ratio`` differ in size, ``theoretical_ratio`` has fewer than two
        different finite values, or ``density_grid`` is not monotonic.
    """
    with np.errstate(divide="ignore", invalid="ignore"):
        quotient = _with_standard_deviation(intensity_numerator).divide(_with_standard_deviation(intensity_denominator))
    ratio = np.where(np.isfinite(quotient.data), quotient.data, np.nan)
    if quotient.mask is not None:
        ratio = np.where(quotient.mask, np.nan, ratio)
    ratio = u.Quantity(ratio, quotient.unit).to(u.one)
    ratio_uncertainty = None
    # The ratios that choose the branch, or the ends of their error bars
    observed = ratio.value[np.newaxis]
    if quotient.uncertainty is not None:
        ratio_uncertainty = u.Quantity(quotient.uncertainty.represent_as(StdDevUncertainty).array, quotient.unit)
        ratio_uncertainty = np.where(np.isnan(ratio), np.nan, ratio_uncertainty.to(u.one))
        observed = u.Quantity([ratio - ratio_uncertainty, ratio + ratio_uncertainty]).value

    density = np.ravel(u.Quantity(density_grid))
    theoretical_ratio = np.ravel(u.Quantity(theoretical_ratio, u.one).value)
    if density.shape != theoretical_ratio.shape:
        msg = "density_grid and theoretical_ratio must have the same size."
        raise ValueError(msg)
    finite = np.isfinite(density.value) & np.isfinite(theoretical_ratio)
    density, theoretical_ratio = density[finite], theoretical_ratio[finite]
    if np.unique(theoretical_ratio).size < 2:
        msg = "theoretical_ratio must have at least two different finite values."
        raise ValueError(msg)
    step = np.diff(density.value)
    if not (np.all(step > 0) or np.all(step < 0)):
        msg = "density_grid must increase or decrease monotonically."
        raise ValueError(msg)

    # The monotonic branches, which share the sample where the curve turns
    branches = []
    start = 0
    previous_sign = 0
    for i, sign in enumerate(np.sign(np.diff(theoretical_ratio))):
        if sign == 0:
            continue
        if previous_sign not in (0, sign):
            branches.append(slice(start, i + 1))
            start = i
        previous_sign = sign
    branches.append(slice(start, None))
    votes = observed[np.isfinite(observed)]

    def score(branch):
        # Ties go to the longer branch, then to the one with the wider range
        curve = theoretical_ratio[branch]
        return ((votes >= curve.min()) & (votes <= curve.max())).sum(), curve.size, np.ptp(curve)

    branch = max(branches, key=score)
    density_grid, theoretical_ratio = density[branch], theoretical_ratio[branch]

    density_lower = density_upper = None
    if ratio_uncertainty is not None:
        # A side of the error bar off the curve gives the curve's end, as in iris_ne_oiv.pro
        curve_ends = density_grid[[0, -1]] if theoretical_ratio[-1] > theoretical_ratio[0] else density_grid[[-1, 0]]
        density_minus, density_plus = _map_ratio_to_quantity(
            observed, density_grid, theoretical_ratio, fill_value=tuple(curve_ends.value)
        )
        density_lower = np.fmin(density_minus, density_plus)
        density_upper = np.fmax(density_minus, density_plus)

    return {
        "ratio": ratio,
        "ratio_uncertainty": ratio_uncertainty,
        "density": _map_ratio_to_quantity(ratio.value, density_grid, theoretical_ratio),
        "density_lower": density_lower,
        "density_upper": density_upper,
        "density_grid": density_grid,
        "theoretical_ratio": u.Quantity(theoretical_ratio, u.one),
    }
