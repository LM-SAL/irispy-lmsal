import astropy.units as u
import numpy as np
import numpy.testing as np_test
import pytest

from irispy.utils.constants import DN_UNIT
from irispy.utils.spectrograph import convert_between_DN_and_photons

# Arrays of DN
SOURCE_DATA_DN = np.array([[0.563, 1.132, -1.343], [-0.719, 1.441, 1.566]])
SOURCE_DATA_DN_1 = np.array([[1, 2, 3], [4, 5, 6]])

# Arrays relating SOURCE_DATA_DN to photons in NUV and FUV
SOURCE_DATA_PHOTONS_NUV = np.array([[10.134, 20.376, -24.174], [-12.942, 25.938, 28.188]])
SOURCE_DATA_PHOTONS_FUV = np.array([[2.252, 4.528, -5.372], [-2.876, 5.764, 6.264]])

# Arrays relating SOURCE_DATA_DN_1 and photons in SJI
SOURCE_DATA_PHOTONS_SJI_1 = np.array([[18, 36, 54], [72, 90, 108]])

single_exposure_time = 2.0
EXPOSURE_TIME = np.zeros(3) + single_exposure_time


@pytest.mark.parametrize(
    "data_arrays, old_unit, new_unit, expected_data_arrays, expected_unit",
    [
        (
            [SOURCE_DATA_DN, SOURCE_DATA_DN],
            DN_UNIT["FUV"],
            u.photon,
            [SOURCE_DATA_PHOTONS_FUV, SOURCE_DATA_PHOTONS_FUV],
            u.photon,
        ),
        (
            [SOURCE_DATA_DN, SOURCE_DATA_DN],
            DN_UNIT["NUV"],
            u.photon,
            [SOURCE_DATA_PHOTONS_NUV, SOURCE_DATA_PHOTONS_NUV],
            u.photon,
        ),
        (
            [SOURCE_DATA_DN_1, SOURCE_DATA_DN_1],
            DN_UNIT["SJI"],
            u.photon,
            [SOURCE_DATA_PHOTONS_SJI_1, SOURCE_DATA_PHOTONS_SJI_1],
            u.photon,
        ),
        (
            [SOURCE_DATA_PHOTONS_FUV, SOURCE_DATA_PHOTONS_FUV],
            u.photon,
            DN_UNIT["FUV"],
            [SOURCE_DATA_DN, SOURCE_DATA_DN],
            DN_UNIT["FUV"],
        ),
        (
            [SOURCE_DATA_PHOTONS_NUV, SOURCE_DATA_PHOTONS_NUV],
            u.photon,
            DN_UNIT["NUV"],
            [SOURCE_DATA_DN, SOURCE_DATA_DN],
            DN_UNIT["NUV"],
        ),
        (
            [SOURCE_DATA_PHOTONS_SJI_1, SOURCE_DATA_PHOTONS_SJI_1],
            u.photon,
            DN_UNIT["SJI"],
            [SOURCE_DATA_DN_1, SOURCE_DATA_DN_1],
            DN_UNIT["SJI"],
        ),
    ],
)
def test_convert_between_DN_and_photons(data_arrays, old_unit, new_unit, expected_data_arrays, expected_unit):
    output_arrays, output_unit = convert_between_DN_and_photons(data_arrays, old_unit, new_unit)
    for i, output_array in enumerate(output_arrays):
        np_test.assert_allclose(output_array, expected_data_arrays[i])
    assert output_unit == expected_unit
