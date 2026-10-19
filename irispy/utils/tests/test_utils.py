import numpy as np
import numpy.testing as np_test
import pytest

from irispy import utils
from irispy.utils.constants import DN_UNIT, READOUT_NOISE


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


@pytest.mark.parametrize(
    ("dtype", "expected"),
    [(np.int16, (4.0, 16.0)), (np.float64, (3.99920004, 9.00360036))],
)
def test_image_clipping(dtype, expected):
    image = np.array([0, 0, 0, 1, 2, 3, 4, 9], dtype=dtype)

    np_test.assert_allclose(utils.image_clipping(image, cutoff=0.2, gamma=2), expected)


def test_import_optional_missing_module_names_the_extra():
    with pytest.raises(ImportError, match=r"pip install 'irispy-lmsal\[rgb\]'"):
        utils.utils._import_optional("definitely_not_a_module", reason="testing", extra="rgb")


def test_import_optional_preserves_transitive_import_error(monkeypatch):
    def import_module(_):
        msg = "No module named 'transitive_dependency'"
        raise ModuleNotFoundError(msg, name="transitive_dependency")

    monkeypatch.setattr(utils.utils, "import_module", import_module)
    with pytest.raises(ModuleNotFoundError, match="transitive_dependency"):
        utils.utils._import_optional("optional_module", reason="testing", extra="tests")


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_calculate_uncertainty_keeps_the_data_dtype(dtype):
    # Level 2 data are float32; the readout noise scalar used to promote the uncertainty to float64
    data = np.array([0.0, 50.0, 3000.0], dtype=dtype)
    uncertainty = utils.calculate_uncertainty(data, READOUT_NOISE["SJI_1400"], DN_UNIT["SJI_1400"])
    assert uncertainty.dtype == dtype
    expected = utils.calculate_uncertainty(data.astype(np.float64), READOUT_NOISE["SJI_1400"], DN_UNIT["SJI_1400"])
    np_test.assert_allclose(uncertainty, expected, rtol=1e-6)


def test_calculate_uncertainty_of_negative_data():
    data = np.array([-200.0, -3.0, 0.0, 50.0, np.nan])
    uncertainty = utils.calculate_uncertainty(data, READOUT_NOISE["FUV"], DN_UNIT["FUV"])
    readout = READOUT_NOISE["FUV"].to(DN_UNIT["FUV"]).value
    np_test.assert_allclose(uncertainty[:3], readout)
    assert uncertainty[3] > readout
    assert np.isnan(uncertainty[4])
