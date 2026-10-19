import numpy as np
import pytest

import astropy.units as u

from irispy.obsid import ObsID

OBSIDS = [
    3400109162,
    3677508065,
    3880903651,
    4050607445,
    4180000056,
]
TEST_DATA = {}
TEST_DATA["exptime"] = [
    8 * u.s,
    8 * u.s,
    30 * u.s,
    4 * u.s,
    1 * u.s,
]
TEST_DATA["raster_desc"] = [
    "Very large coarse 64-step raster",
    "Very large dense 96-step raster",
    "Small sit-and-stare",
    "Very large dense raster (tight timing)",
    "Large coarse raster (tight timing)",
]
TEST_DATA["sjis"] = [
    "C II   Si IV   Mg II h/k",
    "C II   Si IV   Mg II h/k   Mg II w",
    "Si IV",
    "Si IV   Mg II h/k",
    "C II   Si IV   Mg II h/k   Mg II w",
]
TEST_DATA["binning"] = [
    "Spatial x 2, Spectral x 2",
    "Spatial x 1, Spectral x 1",
    "Spatial x 2, Spectral x 8",
    "Spatial x 2, Spectral x 2",
    "Spatial x 1, Spectral x 1",
]
TEST_DATA["fuv_binning"] = [
    "FUV binned same as NUV",
    "FUV spectrally rebinned x 4",
    "FUV spectrally rebinned x 8",
    "FUV spectrally rebinned x 4",
    "FUV binned same as NUV",
]
TEST_DATA["sji_cadence"] = [
    "SJI cadence default",
    "SJI cadence 0.5x faster",
    "SJI cadence default",
    "SJI cadence 10s",
    "SJI cadence 10s",
]
TEST_DATA["linelist"] = [
    "Large linelist",
    "Flare linelist 1",
    "Full readout",
    "Small linelist",
    "Full readout",
]


@pytest.mark.parametrize(
    ("attr_name", "test_input", "expected_output"),
    [(name, obs, output[i]) for (name, output) in TEST_DATA.items() for i, obs in enumerate(OBSIDS)],
)
def test_attribute(attr_name, test_input, expected_output):
    assert ObsID(test_input)[attr_name] == expected_output


@pytest.mark.parametrize("obsid", OBSIDS)
def test_options_keep_description_whitespace(obsid):
    decoded = ObsID(obsid)
    options = decoded.options
    assert options["sjis"]["C II   Si IV   Mg II h/k   Mg II w   "] == 0
    assert options["exptime"][1 * u.s] == 0
    assert type(decoded["raster_fov"]) is str
    assert type(decoded["raster_step"]) is (str if str(obsid).startswith(("40", "41")) else np.float64)
    assert all(type(key) is str and type(value) is int for key, value in options["sjis"].items())


@pytest.mark.parametrize(
    ("test_input", "message"),
    [
        (4204700138, "calibration or special sequence"),
        (4100000010, "no combination of v40 OBS table options"),
        (4050607495, "last two numbers must be between 10 and 72"),
        (3880903650, "last two numbers must be between 1 and 99"),
        (3680903685, "last two numbers must be between 1 and 80"),
        (335987081297, "must have 10 digits"),
        (40, "must have 10 digits"),
    ],
)
def test_invalid_obsid(test_input, message):
    with pytest.raises(ValueError, match=message):
        ObsID(test_input)
