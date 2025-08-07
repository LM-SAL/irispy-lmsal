import numpy as np
import numpy.testing as np_test
import pytest
from matplotlib import pyplot as plt

from astropy import units as u

from sunpy.time import parse_time

from irispy import utils
from irispy.tests.helpers import figure_test

data_dust = np.array(
    [
        [[-1, 2, -3, 4], [2, -200, 5, 3], [0, 1, 2, -300]],
        [[2, -200, 5, 1], [10, -5, 2, 2], [10, -3, 3, 0]],
    ],
)
dust_mask_expected = np.array(
    [
        [
            [True, True, True, True],
            [True, True, True, True],
            [True, True, False, False],
        ],
        [[True, True, True, False], [True, True, True, True], [True, True, True, True]],
    ],
)


@pytest.mark.parametrize(
    ("test_input", "expected_output"),
    [
        ({"detector type": "FUV1"}, "FUV"),
        ({"detector type": "NUV"}, "NUV"),
        ({"detector type": "SJI"}, "SJI"),
    ],
)
def test_get_detector_type(test_input, expected_output):
    assert utils.get_detector_type(test_input) == expected_output


@pytest.mark.parametrize(("input_array", "expected_array"), [(data_dust, dust_mask_expected)])
def test_calculate_dust_mask(input_array, expected_array):
    np_test.assert_array_equal(utils.calculate_dust_mask(input_array), expected_array)


@figure_test
def test_plot_get_interpolated_effective_area():
    # No idea if this is correct as of yet.
    start_obs = parse_time("2013-07-20T17:10:23")
    obs_wavelength = np.linspace(1400.5, 1404.9915000926703, num=692, endpoint=True)
    effective_area = utils.get_interpolated_effective_area(
        start_obs,
        detector_type="FUV",
        obs_wavelength=obs_wavelength * u.Angstrom,
    )
    assert effective_area.shape == obs_wavelength.shape

    fig, ax = plt.subplots()
    ax.plot(obs_wavelength, effective_area)
    ax.set_xlabel("Wavelength (Angstrom)")
    ax.set_ylabel("Effective Area (AA^2) ?")
    ax.set_title("Interpolated Effective Area")
    ax.grid()
    return fig
