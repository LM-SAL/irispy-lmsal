import numpy as np
import pytest

import astropy.units as u

from irispy.utils.lines import get_lines


@pytest.mark.parametrize("unit", [u.angstrom, u.nm, u.Hz, u.eV])
def test_wavelength_range(unit):
    table = get_lines()
    low, high = table["wavelength"][[100, 200]]
    # Reversed in Å, in order once converted to frequency or energy; both ends are included.
    result = get_lines(u.Quantity([high, low]).to(unit, equivalencies=u.spectral()))
    assert result["wavelength"][0] == low
    assert result["wavelength"][-1] == high
    assert len(result) == np.count_nonzero((table["wavelength"] >= low) & (table["wavelength"] <= high))
    assert np.all(np.diff(result["wavelength"]) >= 0)


def test_ranking():
    table = get_lines(region="flare", abundance="coronal", minimum_strength=0)
    assert len(table) == len(get_lines())
    band_order = np.array([("FUV1", "FUV2", "NUV").index(band) for band in table["passband"]])
    assert np.all(np.diff(band_order) >= 0)
    for band in ("FUV1", "FUV2", "NUV"):
        values = table["intensity_flare_coronal"][table["passband"] == band]
        ranked = np.isfinite(values)
        # Strongest first, then the lines without a prediction.
        assert np.all(np.diff(values[ranked]) <= 0)
        assert np.all(ranked[: ranked.sum()])
        assert not np.any(ranked[ranked.sum() :])


def test_minimum_strength():
    everything = get_lines(region="active_region", abundance="photospheric", minimum_strength=0)
    strength = everything["intensity_active_region_photospheric"]
    keep = everything["main"] | ((everything["ion_stage"] >= 3) & (strength >= 0.02))
    # Undocumented cool ions such as Cr II pass the threshold but are left out.
    assert np.any(~keep & (strength >= 0.02))
    short = get_lines(region="active_region", abundance="photospheric")
    np.testing.assert_array_equal(short["wavelength"], everything["wavelength"][keep])


def test_main_only():
    documented = get_lines(main_only=True)
    assert len(documented) == 25
    assert np.all(documented["reference"] != "")


def test_uncategorized_lines():
    # Uncategorized lines hold empty strings, not masked values.
    assert set(get_lines()["category"]) == {"", "flare", "coronal", "transition_region", "chromospheric", "cool_metal"}


def test_results_are_independent():
    changed = get_lines()
    changed["wavelength"][0] = 1 * u.angstrom
    changed.meta["passband_limits_angstrom"]["FUV1"][0] = 0
    fresh = get_lines()
    assert fresh["wavelength"][0] != 1 * u.angstrom
    assert fresh.meta["passband_limits_angstrom"]["FUV1"][0] != 0


@pytest.mark.parametrize(
    "kwargs",
    [
        {"region": "sunspot", "abundance": "coronal"},
        {"region": "quiet_sun"},
        {"abundance": "coronal"},
        {"abundance": "feldman"},
        {"minimum_strength": 0.1},
        {"wavelength_range": [1, 2, 3] * u.angstrom},
        {"wavelength_range": [0, 1] * u.angstrom},
        {"wavelength_range": [np.nan, 1] * u.angstrom},
        {"wavelength_range": [np.inf, 1] * u.eV},
        {"wavelength_range": 1400 * u.angstrom},
    ],
)
def test_invalid_arguments(kwargs):
    with pytest.raises(ValueError, match=r"region|abundance|minimum_strength|wavelength_range"):
        get_lines(**kwargs)


@pytest.mark.parametrize("wavelength_range", [[1, 2] * u.s, [1402, 1404]])
def test_non_spectral_units(wavelength_range):
    with pytest.raises(u.UnitConversionError, match="wavelength_range must be in a spectral unit"):
        get_lines(wavelength_range)
