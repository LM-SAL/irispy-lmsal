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
    4191006128,  # 20130718_213014, a sparse/dense/sparse raster; matches IDL iris_obsid2struct
]
TEST_DATA = {}
TEST_DATA["exptime"] = [
    8 * u.s,
    8 * u.s,
    30 * u.s,
    4 * u.s,
    1 * u.s,
    4 * u.s,
]
TEST_DATA["raster_desc"] = [
    "Very large coarse 64-step raster",
    "Very large dense 96-step raster",
    "Small sit-and-stare",
    "Very large dense raster (tight timing)",
    "Large coarse raster (tight timing)",
    "Large sparse/dense/sparse raster",
]
TEST_DATA["sjis"] = [
    "C II   Si IV   Mg II h/k",
    "C II   Si IV   Mg II h/k   Mg II w",
    "Si IV",
    "Si IV   Mg II h/k",
    "C II   Si IV   Mg II h/k   Mg II w",
    "C II   Si IV   Mg II h/k   Mg II w s",
]
TEST_DATA["binning"] = [
    "Spatial x 2, Spectral x 2",
    "Spatial x 1, Spectral x 1",
    "Spatial x 2, Spectral x 8",
    "Spatial x 2, Spectral x 2",
    "Spatial x 1, Spectral x 1",
    "Spatial x 1, Spectral x 1",
]
TEST_DATA["fuv_binning"] = [
    "FUV binned same as NUV",
    "FUV spectrally rebinned x 4",
    "FUV spectrally rebinned x 8",
    "FUV spectrally rebinned x 4",
    "FUV binned same as NUV",
    "FUV binned same as NUV",
]
TEST_DATA["sji_cadence"] = [
    "SJI cadence default",
    "SJI cadence 0.5x faster",
    "SJI cadence default",
    "SJI cadence 10s",
    "SJI cadence 10s",
    "SJI cadence 0.25x faster",
]
TEST_DATA["linelist"] = [
    "Large linelist",
    "Flare linelist 1",
    "Full readout",
    "Small linelist",
    "Full readout",
    "Full readout",
]
TEST_DATA["raster_step"] = [
    2 * u.arcsec,
    0.33 * u.arcsec,
    0.33 * u.arcsec,
    0.33 * u.arcsec,
    2 * u.arcsec,
    [1, 0.33, 1] * u.arcsec,
]
TEST_DATA["spec_cadence"] = [
    64 * u.s,
    96 * u.s,
    1 * u.s,
    400 * u.s,
    64 * u.s,
    60 * u.s,
]
TEST_DATA["raster_fulldesc"] = [
    "Very large coarse 64-step raster 126x175 64s",
    "Very large dense 96-step raster 31.35x175 96s",
    "Small sit-and-stare 0.3x30 1s",
    "Very large dense raster (tight timing) 132x175 400s",
    "Large coarse raster (tight timing) 126x120 64s",
    "Large sparse/dense/sparse raster 11x120  11.67x120  11x120  t: 33.67x120 60s",
]
TEST_DATA["sji_fov"] = ["175x175", "175x175", "60x60", "175x175", "120x120", "120x120"]


@pytest.mark.parametrize(
    ("attr_name", "test_input", "expected_output"),
    [(name, obs, output[i]) for (name, output) in TEST_DATA.items() for i, obs in enumerate(OBSIDS)],
)
def test_attribute(attr_name, test_input, expected_output):
    assert np.all(ObsID(test_input)[attr_name] == expected_output)


@pytest.mark.parametrize("obsid", OBSIDS)
def test_values_have_one_type(obsid):
    # FITS headers store OBSID as a string
    decoded = ObsID(str(obsid))
    assert type(decoded["obsid"]) is int
    for key in decoded.keys() - {"obsid", "raster_step", "spec_cadence", "exptime"}:
        assert type(decoded[key]) is str
    assert decoded["raster_step"].shape == ((3,) if obsid == 4191006128 else ())
    assert decoded["spec_cadence"].unit == u.s
    for key, options in decoded.options.items():
        assert decoded[key] in options


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
