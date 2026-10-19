import types

import numpy as np
import pooch
import pytest

import astropy.units as u

from irispy.utils.density import density_diagnostic, map_ratio_to_quantity


def _fake_ion(ratio, **attributes):
    # Contribution functions of two transitions in the ratio ``ratio``
    calls = {}

    def contribution_function(density, **kwargs):
        calls["density"] = density
        calls["kwargs"] = kwargs
        return np.stack([ratio, np.ones_like(ratio)], axis=-1) * u.erg * u.cm**3 / u.s

    transitions = types.SimpleNamespace(wavelength=[1399.776, 1401.163] * u.AA, is_bound_bound=np.array([True, True]))
    return types.SimpleNamespace(
        transitions=transitions, contribution_function=contribution_function, calls=calls, **attributes
    )


def test_map_ratio_to_quantity():
    density = [1e8, 1e9, 1e10] * u.cm**-3
    theoretical_ratio = [0.2, 0.5, 0.8] * u.dimensionless_unscaled

    mapped_density = map_ratio_to_quantity([0.2, 0.65, 1.0], density, theoretical_ratio)

    assert mapped_density.unit == density.unit
    assert mapped_density.shape == (3,)
    assert u.allclose(mapped_density[:2], [1e8, 5.5e9] * u.cm**-3)
    assert np.isnan(mapped_density[-1].value)


def test_map_ratio_to_quantity_decreasing_curve():
    density = [1e8, 1e9, 1e10] * u.cm**-3
    theoretical_ratio = [0.8, 0.5, 0.2] * u.dimensionless_unscaled

    assert u.allclose(map_ratio_to_quantity(0.35, density, theoretical_ratio), 5.5e9 * u.cm**-3)


@pytest.mark.parametrize(
    ("observed_ratio", "theoretical_ratio", "kwargs", "message"),
    [
        (0.35, [0.2, 0.8], {}, "must have the same shape"),
        (0.35, [0.2, np.nan, np.nan], {}, "At least two finite samples"),
        (0.35, [0.5, 0.5, 0.5], {}, "must vary over the provided quantity grid"),
        (0.35, [0.2, 0.8, 0.5], {}, "must be monotonic"),
        (1.5, [0.2, 0.5, 0.8], {"bounds_error": True}, "outside the theoretical ratio range"),
    ],
    ids=["shape", "finite-samples", "constant", "non-monotonic", "bounds"],
)
def test_map_ratio_to_quantity_rejects_bad_curves(observed_ratio, theoretical_ratio, kwargs, message):
    density = [1e8, 1e9, 1e10] * u.cm**-3

    with pytest.raises(ValueError, match=message):
        map_ratio_to_quantity(observed_ratio, density, theoretical_ratio * u.dimensionless_unscaled, **kwargs)


@pytest.mark.parametrize("shape", [(), (2,), (1, 2)])
def test_density_diagnostic(shape):
    fake_ion = _fake_ion([0.5, 1.0], temperature=[1e5, 2e5] * u.K, formation_temperature=2e5 * u.K)

    result = density_diagnostic(
        np.resize([1.2, 1.5], shape) * u.ct,
        np.full(shape, 2.0) * u.ct,
        [1e10, 1e11] * u.cm**-3,
        ion=fake_ion,
        numerator=1399.78 * u.angstrom,
        denominator=1401.16 * u.angstrom,
        intensity_numerator_uncertainty=np.full(shape, 0.1) * u.ct,
        intensity_denominator_uncertainty=np.full(shape, 0.2) * u.ct,
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


def test_density_diagnostic_builds_theoretical_ratio_from_the_ion():
    fake_ion = _fake_ion(
        [
            [0.1, 0.2, 0.3],
            [0.2, 0.4, 0.6],
            [0.3, 0.6, 0.9],
        ],
        temperature=[1e5, 2e5, 3e5] * u.K,
        formation_temperature=2e5 * u.K,
    )
    density_grid = [1, 2, 3] * u.cm**-3

    result = density_diagnostic(
        [0.4],
        [1.0],
        density_grid,
        ion=fake_ion,
        numerator=1399.78 * u.angstrom,
        denominator=1401.16 * u.angstrom,
        line_ratio_kwargs={"use_two_ion_model": False},
    )

    assert u.allclose(fake_ion.calls["density"], density_grid)
    assert fake_ion.calls["kwargs"] == {"use_two_ion_model": False}
    np.testing.assert_allclose(result["theoretical_ratio"].value, [0.2, 0.4, 0.6])
    assert u.allclose(result["density"], [2] * u.cm**-3)


def test_density_diagnostic_error_bar_off_the_curve_gives_its_end():
    # The O IV 1399.78 / 1401.16 Å curve of SolarSoft's iris_ne_oiv.pro; 0.375 + 0.078 lies above it
    oiv_ratio = [
        0.16710843, 0.1670735, 0.16709631, 0.16718450, 0.16736760, 0.16770638, 0.16831132, 0.16937367, 0.17121214,
        0.17433385, 0.17949909, 0.18777250, 0.20053728, 0.2193597, 0.24536467, 0.2779118, 0.3135369, 0.34697879,
        0.3740375, 0.39339120, 0.40603911, 0.41382271, 0.4184366, 0.42111109, 0.42264118,
    ]  # fmt: skip

    result = density_diagnostic(
        0.375,
        1.0,
        np.logspace(7, 13, 25) * u.cm**-3,
        ion=_fake_ion(oiv_ratio),
        numerator=1399.78 * u.angstrom,
        denominator=1401.16 * u.angstrom,
        intensity_numerator_uncertainty=0.078,
        intensity_denominator_uncertainty=0.0,
    )

    assert result["density_upper"] == 1e13 * u.cm**-3
    # iris_ne_oiv.pro gives log densities 11.523 and 10.889 from the nearest of 200 points, 0.03 dex apart,
    # on a spline of its curve
    np.testing.assert_allclose(
        np.log10(u.Quantity([result["density"], result["density_lower"]]).to_value(u.cm**-3)),
        [11.523, 10.889],
        atol=0.02,
    )


def test_density_diagnostic_adds_transitions_within_0_02_angstrom():
    # CHIANTI's O IV 1401.149 and 1401.163 Å transitions both lie within 0.02 Å of the catalog's 1401.157 Å
    ion = types.SimpleNamespace(
        transitions=types.SimpleNamespace(
            wavelength=[1399.776, 1401.163, 1401.149] * u.AA, is_bound_bound=np.array([True, True, True])
        ),
        contribution_function=lambda _density: [[1.0, 0.3, 0.1], [1.0, 0.6, 0.2]] * u.erg * u.cm**3 / u.s,
    )

    result = density_diagnostic(
        0.6,
        1.0,
        [1e10, 1e11] * u.cm**-3,
        ion=ion,
        numerator=1401.157 * u.angstrom,
        denominator=1399.776 * u.angstrom,
    )

    np.testing.assert_allclose(result["theoretical_ratio"].value, [0.3 + 0.1, 0.6 + 0.2])


def test_density_diagnostic_rejects_wavelength_without_transition():
    with pytest.raises(ValueError, match=r"No transition .* within 0\.02 Å .* nearest are at 1399\.776, 1401\.163 Å"):
        density_diagnostic(
            1.0,
            1.0,
            [1e10, 1e11] * u.cm**-3,
            ion=_fake_ion([0.5, 1.0]),
            numerator=1404.806 * u.angstrom,
            denominator=1401.163 * u.angstrom,
        )


@pytest.mark.remote_data
def test_density_diagnostic_oiv_curve_from_fiasco():
    fiasco = pytest.importorskip("fiasco")
    # The CHIANTI database of the documentation, with O IV only
    database = pooch.retrieve(
        "https://github.com/LM-SAL/irispy-data/releases/download/v1/chianti_dbase_o4.h5",
        known_hash="7a7a6eb1bfd043695b6f52c641b663d71f8a33a836da35469f54151af57d014d",
    )
    ion = fiasco.Ion("O IV", 10**5.15 * u.K, hdf5_dbase_root=database)

    result = density_diagnostic(
        1.0,
        1.0,
        np.logspace(9, 12, 7) * u.cm**-3,
        ion=ion,
        numerator=1399.776 * u.angstrom,
        denominator=1401.157 * u.angstrom,
        line_ratio_kwargs={"use_two_ion_model": False},
    )

    # fiasco.line_ratio of fiasco 0.9.dev4 with this database and the 1401.149 and 1401.163 Å transitions;
    # fiasco 0.8.2 gives up to 0.006% less
    expected = [0.1720960770, 0.1803754352, 0.2014199870, 0.2462960737, 0.3144225710, 0.3746827316, 0.4064668318]
    np.testing.assert_allclose(result["theoretical_ratio"].value, expected, rtol=1e-4)


def test_density_diagnostic_requires_ion():
    with pytest.raises(ValueError, match="ion, numerator, and denominator are required"):
        density_diagnostic(
            [4, 8] * u.ct,
            [2, 4] * u.ct,
            [1e10, 1e11] * u.cm**-3,
            ion=None,
            numerator=1399.78 * u.angstrom,
            denominator=1401.16 * u.angstrom,
        )


def test_density_diagnostic_selects_monotonic_branch():
    fake_ion = _fake_ion([0.1, 0.2, 0.3, 0.25, 0.2], temperature=[1e5, 2e5] * u.K, formation_temperature=2e5 * u.K)

    result = density_diagnostic(
        [2.4, 2.8],
        [10.0, 10.0],
        density_grid=[1, 2, 3, 4, 5] * u.cm**-3,
        ion=fake_ion,
        numerator=1399.78 * u.angstrom,
        denominator=1401.16 * u.angstrom,
    )

    np.testing.assert_allclose(result["theoretical_ratio"].value, [0.1, 0.2, 0.3])
    np.testing.assert_allclose(result["density_grid"].value, [1, 2, 3])


def test_density_diagnostic_zero_denominator():
    fake_ion = _fake_ion([0.5, 1.0], temperature=[1e5, 2e5] * u.K, formation_temperature=2e5 * u.K)

    result = density_diagnostic(
        [1.0, 1.5] * u.ct,
        [2.0, 0.0] * u.ct,
        [1e10, 1e11] * u.cm**-3,
        ion=fake_ion,
        numerator=1399.78 * u.angstrom,
        denominator=1401.16 * u.angstrom,
    )

    np.testing.assert_allclose(result["ratio"].value, [0.5, np.nan])
    assert result["ratio_uncertainty"] is None
    assert result["density_lower"] is None
    assert result["density_upper"] is None


def test_density_diagnostic_partial_uncertainty_raises():
    with pytest.raises(ValueError, match="Both numerator and denominator uncertainties"):
        density_diagnostic(
            [4, 8] * u.ct,
            [2, 4] * u.ct,
            [1e10, 1e11] * u.cm**-3,
            ion=types.SimpleNamespace(),
            numerator=1399.78 * u.angstrom,
            denominator=1401.16 * u.angstrom,
            intensity_numerator_uncertainty=[0.1, 0.1] * u.ct,
        )


def test_density_diagnostic_explicit_temperature():
    """
    Regression: explicit scalar temperature must not be clobbered by
    formation_temperature.
    """
    fake_ion = _fake_ion(
        [
            [0.1, 0.2, 0.3],
            [0.2, 0.4, 0.6],
            [0.3, 0.6, 0.9],
        ],
        temperature=[1e5, 2e5, 3e5] * u.K,
        formation_temperature=1e5 * u.K,
    )
    density_grid = [1, 2, 3] * u.cm**-3

    result = density_diagnostic(
        [0.4],
        [1.0],
        density_grid,
        ion=fake_ion,
        numerator=1399.78 * u.angstrom,
        denominator=1401.16 * u.angstrom,
        temperature=2e5 * u.K,
    )

    # If temperature=2e5 K were ignored, the fallback formation temperature
    # would use the first row and return density=1.
    assert u.allclose(result["density"], [2] * u.cm**-3)


def test_density_diagnostic_rejects_vector_temperature():
    fake_ion = _fake_ion(
        [
            [0.1, 0.2, 0.3],
            [0.2, 0.4, 0.6],
            [0.3, 0.6, 0.9],
        ],
        temperature=[1e5, 2e5, 3e5] * u.K,
        formation_temperature=1e5 * u.K,
    )

    with pytest.raises(ValueError, match="temperature must be scalar"):
        density_diagnostic(
            [0.4],
            [1.0],
            [1, 2, 3] * u.cm**-3,
            ion=fake_ion,
            numerator=1399.78 * u.angstrom,
            denominator=1401.16 * u.angstrom,
            temperature=[1e5, 2e5] * u.K,
        )


def test_density_diagnostic_uncertainty_shape_mismatch():
    fake_ion = _fake_ion([0.5, 1.0], temperature=[1e5, 2e5] * u.K, formation_temperature=2e5 * u.K)

    with pytest.raises(ValueError, match="Uncertainty shapes"):
        density_diagnostic(
            [1.0, 1.5] * u.ct,
            [2.0, 2.0] * u.ct,
            [1e10, 1e11] * u.cm**-3,
            ion=fake_ion,
            numerator=1399.78 * u.angstrom,
            denominator=1401.16 * u.angstrom,
            intensity_numerator_uncertainty=[[0.1], [0.1]] * u.ct,
            intensity_denominator_uncertainty=[[0.2], [0.2]] * u.ct,
        )


def test_density_diagnostic_rejects_ratio_with_other_temperatures():
    # Three temperature rows from fiasco for an ion with two temperatures
    fake_ion = _fake_ion(
        [[0.1, 0.2], [0.2, 0.4], [0.3, 0.6]], temperature=[1e5, 2e5] * u.K, formation_temperature=2e5 * u.K
    )

    with pytest.raises(ValueError, match=r"temperature axis does not match ion\.temperature"):
        density_diagnostic(
            [0.4],
            [1.0],
            [1, 2] * u.cm**-3,
            ion=fake_ion,
            numerator=1399.78 * u.angstrom,
            denominator=1401.16 * u.angstrom,
        )
