"""
``cut.py`` keeps the world coordinates and the values of the pixels it keeps.
"""

import numpy as np
from cut import cut  # tools/ is on sys.path while pytest collects this file

from astropy.io import fits

from irispy.data.test import get_test_filepath

SJI = "sns/iris_l2_20210905_001833_3620258102_SJI_2832_t000_test.fits"
RASTER = (
    "raster/iris_l2_20140329_140938_3860258481_raster/iris_l2_20140329_140938_3860258481_raster_t000_r00000_test.fits"
)


def _intermediate(header, axis, pixels):
    """
    The linear coordinate of 0-based ``pixels`` along ``axis``, before the PC matrix and
    the projection.
    """
    return (np.asarray(pixels) + 1 - header[f"CRPIX{axis}"]) * header[f"CDELT{axis}"]


def _assert_same_world(cut_header, kept, source_header, source_pixels):
    """
    The ``kept`` pixels of the cut HDU have the world coordinates of ``source_pixels``
    in the source.
    """
    for axis, (new, old) in enumerate(zip(kept, source_pixels, strict=True), start=1):
        np.testing.assert_allclose(_intermediate(cut_header, axis, new), _intermediate(source_header, axis, old))
    for key in ("CRVAL1", "CRVAL2", "CRVAL3", "PC1_1", "PC1_2", "PC2_1", "PC2_2"):
        assert cut_header[key] == source_header[key]


def test_cut_sji_keeps_world_coordinates_and_values(tmp_path):
    source = get_test_filepath(SJI)
    output = tmp_path / "sji.fits"
    exposures = np.array([1, 4, 7])  # evenly spaced, so the time axis is rescaled too
    rows, columns = slice(5, 30), slice(3, 20)

    cut(source, output, exposures=exposures, rows=rows, columns=columns)

    with fits.open(source, do_not_scale_image_data=True) as src, fits.open(output, do_not_scale_image_data=True) as out:
        np.testing.assert_array_equal(out[0].data, src[0].data[exposures][:, rows, columns])
        assert out[0].header["NEXP"] == 3
        kept = (np.arange(17), np.arange(25), np.arange(3))
        _assert_same_world(out[0].header, kept, src[0].header, (np.arange(3, 20), np.arange(5, 30), exposures))
        np.testing.assert_array_equal(out[-2].data, src[-2].data[exposures])
        assert out[-1].data.tolist() == src[-1].data[exposures].tolist()  # the level 1 filename table


def test_cut_raster_keeps_world_coordinates_and_values(tmp_path):
    source = get_test_filepath(RASTER)
    output = tmp_path / "raster.fits"
    exposures = np.array([2, 3, 4])
    rows, pixels = slice(10, 50), slice(5, 20)
    windows = [("Si IV 1403", pixels, None), ("Mg II k 2796", None, ("Mg II h 2803", 2803.53))]

    cut(source, output, exposures=exposures, rows=rows, windows=windows)

    with fits.open(source, do_not_scale_image_data=True) as src, fits.open(output, do_not_scale_image_data=True) as out:
        primary = out[0].header
        assert primary["NWIN"] == 2
        assert primary["NEXP"] == primary["NRASTERP"] == 3
        assert (primary["TDESC1"], primary["TDESC2"], primary["TWAVE2"]) == ("Si IV 1403", "Mg II h 2803", 2803.53)
        # Si IV 1403 is window 5 of the source and Mg II k 2796 window 9
        np.testing.assert_array_equal(out[1].data, src[5].data[exposures][:, rows, pixels])
        np.testing.assert_array_equal(out[2].data, src[9].data[exposures][:, rows])
        kept = (np.arange(15), np.arange(40), np.arange(3))
        _assert_same_world(out[1].header, kept, src[5].header, (np.arange(5, 20), np.arange(10, 50), exposures))
        kept = (np.arange(63), np.arange(40), np.arange(3))
        _assert_same_world(out[2].header, kept, src[9].header, (np.arange(63), np.arange(10, 50), exposures))
        wavelength = _intermediate(out[1].header, 1, [0, 14]) + out[1].header["CRVAL1"]
        np.testing.assert_allclose([primary["TWMIN1"], primary["TWMAX1"]], wavelength)
        np.testing.assert_array_equal(out[-2].data, src[-2].data[exposures])
        # cut.py repairs the understated TFIELDS of the filename table, so the source rows parse differently
        assert out[-1].header["TFIELDS"] == sum(1 for key in out[-1].header if key.startswith("TTYPE"))
        assert [row[0] for row in out[-1].data] == [row[0] for row in src[-1].data[exposures]]
