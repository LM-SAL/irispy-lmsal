import numpy as np
import numpy.testing as np_test
import pytest

from irispy import utils

DETECTOR_TYPE_KEY = "detector type"


# Arrays for the calcuate_dust_mask method.
data_dust = np.array(
    [
        [[-1, 2, -3, 4], [2, -200, 5, 3], [0, 1, 2, -300]],
        [[2, -200, 5, 1], [10, -5, 2, 2], [10, -3, 3, 0]],
    ]
)
dust_mask_expected = np.array(
    [
        [
            [True, True, True, True],
            [True, True, True, True],
            [True, True, False, False],
        ],
        [[True, True, True, False], [True, True, True, True], [True, True, True, True]],
    ]
)


@pytest.mark.parametrize(
    "test_input, expected_output",
    [
        ({DETECTOR_TYPE_KEY: "FUV1"}, "FUV"),
        ({DETECTOR_TYPE_KEY: "NUV"}, "NUV"),
        ({DETECTOR_TYPE_KEY: "SJI"}, "SJI"),
    ],
)
def test_get_detector_type(test_input, expected_output):
    assert utils.get_detector_type(test_input) == expected_output


@pytest.mark.parametrize("input_array, expected_array", [(data_dust, dust_mask_expected)])
def test_calculate_dust_mask(input_array, expected_array):
    np_test.assert_array_equal(utils.calculate_dust_mask(input_array), expected_array)
