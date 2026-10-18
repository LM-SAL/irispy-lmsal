import re

import pytest

from astropy.io import fits

from irispy.io import read_files, read_sji_lvl2
from irispy.io._level2 import _classify
from irispy.io.spectrograph import read_spectrograph_lvl2
from irispy.io.tests.test_mosaic import HEADER as MOSAIC_HEADER

# The cards of a level 3 file made by iris_make_fits_level3: INSTRUME stays "IRIS", the level is 3.
LEVEL3_CARDS = {"INSTRUME": "IRIS", "DATA_LEV": 3.0, "LVL_NUM": 3.0, "OBSID": "3630104144"}


@pytest.fixture
def level3_file(tmp_path):
    filename = tmp_path / "iris_l3_20150917_073915_3630104144_t000_HI6563_im.fits"
    fits.PrimaryHDU(header=fits.Header(LEVEL3_CARDS)).writeto(filename)
    return filename


@pytest.mark.parametrize(
    ("cards", "expected"),
    [
        ({"INSTRUME": "SJI", "DATA_LEV": 2.0, "LVL_NUM": 2.0}, ("SJI", 2.0)),
        ({"INSTRUME": "AIA_3", "DATA_LEV": 2.0, "LVL_NUM": 2.0}, ("SJI", 2.0)),
        ({"INSTRUME": "SOT-SP", "LVL_NUM": 2.0}, ("SJI", 2.0)),
        ({"INSTRUME": "SPEC", "DATA_LEV": 2.0, "LVL_NUM": 2.0}, ("raster", 2.0)),
        (LEVEL3_CARDS, ("SJI", 3.0)),
        (MOSAIC_HEADER, (None, None)),
        ({"TELESCOP": "SDO/AIA", "INSTRUME": "AIA_3", "LVL_NUM": 1.0}, ("SJI", 1.0)),
        ({}, (None, None)),
    ],
    ids=["sji", "aia", "sot", "raster", "level3", "mosaic", "aia-level1", "empty"],
)
def test_classify(cards, expected):
    assert _classify(fits.Header(cards)) == expected


@pytest.mark.parametrize("kind", ["directory", "text file", "level 3 file"])
def test_read_files_rejects_what_is_not_level_2(tmp_path, level3_file, kind):
    path = {"directory": tmp_path, "text file": tmp_path / "notes.txt", "level 3 file": level3_file}[kind]
    if kind == "text file":
        path.write_text("not a FITS file")

    with pytest.raises(ValueError, match=re.escape(f"{path} is not an IRIS level 2 file")):
        read_files(path)
    with pytest.raises(ValueError, match="No supported IRIS files were loaded"):
        read_files(path, allow_errors=True)


@pytest.mark.parametrize(
    ("reader", "product", "wrong_fixture"),
    [(read_sji_lvl2, "SJI", "sns_sg_file"), (read_spectrograph_lvl2, "raster", "sns_sji_1330_file")],
    ids=["sji-reader", "spectrograph-reader"],
)
def test_readers_reject_other_products_and_levels(request, level3_file, reader, product, wrong_fixture):
    wrong_file = request.getfixturevalue(wrong_fixture)
    for path in (wrong_file, level3_file):
        with pytest.raises(ValueError, match=re.escape(f"{path} is not an IRIS level 2 {product} file")):
            reader(path)


def test_read_spectrograph_lvl2_names_the_wrong_file_in_a_list(sns_sg_file, sns_sji_1330_file):
    with pytest.raises(ValueError, match=re.escape(f"{sns_sji_1330_file} is not an IRIS level 2 raster file")):
        read_spectrograph_lvl2([sns_sg_file, sns_sji_1330_file])
