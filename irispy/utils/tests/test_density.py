import numpy as np
import pytest

import astropy.units as u
from astropy.nddata import NDData, StdDevUncertainty
from astropy.wcs import WCS

from ndcube import NDCube

from irispy.utils.density import density_diagnostic


@pytest.mark.parametrize(
    "theoretical_ratio", [[0.2, 0.5, 0.8], [[0.2, 0.5, 0.8]]], ids=["one-dimensional", "one-temperature-row"]
)
def test_density_diagnostic_interpolates_the_curve(theoretical_ratio):
    density = [1e8, 1e9, 1e10] * u.cm**-3

    result = density_diagnostic([0.2, 0.65, 1.0], 1.0, density, theoretical_ratio)

    assert result["density"].unit == density.unit
    assert u.allclose(result["density"][:2], [1e8, 5.5e9] * u.cm**-3)
    assert np.isnan(result["density"][-1].value)


def test_density_diagnostic_decreasing_curve():
    result = density_diagnostic(0.35, 1.0, [1e8, 1e9, 1e10] * u.cm**-3, [0.8, 0.5, 0.2])

    assert u.allclose(result["density"], 5.5e9 * u.cm**-3)


@pytest.mark.parametrize(
    ("theoretical_ratio", "message"),
    [
        ([0.2, 0.8], "must have the same size"),
        ([0.2, np.nan, np.nan], "at least two different finite values"),
        ([0.5, 0.5, 0.5], "at least two different finite values"),
    ],
    ids=["size", "finite-samples", "constant"],
)
def test_density_diagnostic_rejects_bad_curves(theoretical_ratio, message):
    with pytest.raises(ValueError, match=message):
        density_diagnostic(0.35, 1.0, [1e8, 1e9, 1e10] * u.cm**-3, theoretical_ratio)


def test_density_diagnostic_rejects_an_unordered_density_grid():
    with pytest.raises(ValueError, match="monotonically"):
        density_diagnostic(0.25, 1.0, [1, 3, 2, 4] * u.cm**-3, [0.1, 0.2, 0.3, 0.4])


@pytest.mark.parametrize("shape", [(), (2,), (1, 2)])
def test_density_diagnostic(shape):
    result = density_diagnostic(
        NDData(np.resize([1.2, 1.5], shape) * u.ct, uncertainty=StdDevUncertainty(np.full(shape, 0.1))),
        NDData(np.full(shape, 2.0) * u.ct, uncertainty=StdDevUncertainty(np.full(shape, 0.2))),
        [1e10, 1e11] * u.cm**-3,
        [0.5, 1.0],
    )

    np.testing.assert_allclose(result["ratio"].value, np.resize([0.6, 0.75], shape))
    # Hand-calculated ratio variances; the linear density curve has slope 1.8e11 cm^-3.
    ratio_error = np.sqrt(np.resize([0.0061, 0.008125], shape)) * u.one
    expected_density = np.resize([2.8e10, 5.5e10], shape) * u.cm**-3
    density_error = ratio_error.value * 1.8e11 * u.cm**-3
    assert u.allclose(result["ratio_uncertainty"], ratio_error)
    assert u.allclose(result["density"], expected_density)
    assert u.allclose(result["density_lower"], expected_density - density_error)
    assert u.allclose(result["density_upper"], expected_density + density_error)
    for key in ("ratio", "ratio_uncertainty", "density", "density_lower", "density_upper"):
        assert result[key].shape == shape


def test_density_diagnostic_takes_uncertainty_and_mask_from_cubes():
    numerator = NDCube(
        [0.0, 1.2, 1.5] * u.ct,
        wcs=WCS(naxis=1),
        uncertainty=StdDevUncertainty([0.1, 0.1, 0.1]),
        mask=[False, False, True],
    )

    result = density_diagnostic(numerator, [2.0, 2.0, 2.0] * u.ct, [1e10, 1e11] * u.cm**-3, [0.0, 1.0])

    np.testing.assert_allclose(result["ratio"].value, [0.0, 0.6, np.nan])
    # A numerator of 0 has a finite error bar, 0.1 / 2
    np.testing.assert_allclose(result["ratio_uncertainty"].value, [0.05, 0.05, np.nan])
    assert u.allclose(result["density_upper"][0], 1.45e10 * u.cm**-3)
    for key in ("density", "density_lower", "density_upper"):
        assert np.isnan(result[key][-1])


def test_density_diagnostic_error_bar_off_the_curve_gives_its_end():
    # The O IV 1399.78 / 1401.16 Å curve of SolarSoft's iris_ne_oiv.pro; 0.375 + 0.078 lies above it
    oiv_ratio = [
        0.16710843, 0.1670735, 0.16709631, 0.16718450, 0.16736760, 0.16770638, 0.16831132, 0.16937367, 0.17121214,
        0.17433385, 0.17949909, 0.18777250, 0.20053728, 0.2193597, 0.24536467, 0.2779118, 0.3135369, 0.34697879,
        0.3740375, 0.39339120, 0.40603911, 0.41382271, 0.4184366, 0.42111109, 0.42264118,
    ]  # fmt: skip

    result = density_diagnostic(
        NDData(0.375, uncertainty=StdDevUncertainty(0.078)), 1.0, np.logspace(7, 13, 25) * u.cm**-3, oiv_ratio
    )

    assert result["density_upper"] == 1e13 * u.cm**-3
    # iris_ne_oiv.pro gives log densities 11.523 and 10.889 from the nearest of 200 points, 0.03 dex apart,
    # on a spline of its curve
    np.testing.assert_allclose(
        np.log10(u.Quantity([result["density"], result["density_lower"]]).to_value(u.cm**-3)),
        [11.523, 10.889],
        atol=0.02,
    )


def test_density_diagnostic_error_bar_off_a_flat_end_gives_the_grid_end():
    result = density_diagnostic(
        NDData(0.25, uncertainty=StdDevUncertainty(0.1)), 1.0, [1, 2, 3, 4] * u.cm**-3, [0.1, 0.2, 0.3, 0.3]
    )

    assert result["density_upper"] == 4 * u.cm**-3


def test_density_diagnostic_selects_monotonic_branch():
    result = density_diagnostic([2.4, 2.8], [10.0, 10.0], [1, 2, 3, 4, 5] * u.cm**-3, [0.1, 0.2, 0.3, 0.25, 0.2])

    np.testing.assert_allclose(result["theoretical_ratio"].value, [0.1, 0.2, 0.3])
    np.testing.assert_allclose(result["density_grid"].value, [1, 2, 3])


def test_density_diagnostic_zero_denominator():
    result = density_diagnostic([1.0, 1.5] * u.ct, [2.0, 0.0] * u.ct, [1e10, 1e11] * u.cm**-3, [0.5, 1.0])

    np.testing.assert_allclose(result["ratio"].value, [0.5, np.nan])
    assert result["ratio_uncertainty"] is None
    assert result["density_lower"] is None
    assert result["density_upper"] is None
