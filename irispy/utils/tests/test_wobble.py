import os

import numpy as np
import pytest
from astropy.io import fits

from irispy.data.test import get_test_filepath
from irispy.utils.wobble import wobble_movie


@pytest.fixture()
def filelist():
    return [
        get_test_filepath("iris_l2_20170502_052551_3893010094_SJI_1330_t000.fits"),
        get_test_filepath("iris_l2_20170502_052551_3893010094_SJI_1400_t000.fits"),
        get_test_filepath("iris_l2_20170502_052551_3893010094_SJI_2796_t000.fits"),
        get_test_filepath("iris_l2_20170502_052551_3893010094_SJI_2832_t000.fits"),
    ]


@pytest.fixture(scope="session")
def fake_long_obs(tmp_path_factory):
    header = fits.getheader(get_test_filepath("iris_l2_20170502_052551_3893010094_SJI_2832_t000.fits"))
    header["STARTOBS"] = "2017-05-02T05:25:51.000"
    header["ENDOBS"] = "2017-05-02T08:25:51.000"
    header["NAXIS3"] = 100
    data = np.random.rand(100, header["NAXIS2"], header["NAXIS1"])
    temp_dir = tmp_path_factory.mktemp("IRIS")
    hdu = fits.PrimaryHDU(data=data, header=header, do_not_scale_image_data=True, scale_back=True)
    fits_file = os.fspath(temp_dir.joinpath("iris_l2_20170502_052551_3893010094_SJI_2832_t000.fits"))
    hdu.writeto(fits_file)
    return [fits_file]


def test_wobble_movie(fake_long_obs, tmp_path):
    movies = wobble_movie(fake_long_obs, outdir=tmp_path)
    assert movies != []
    movies = wobble_movie(fake_long_obs, outdir=tmp_path, trim=True)
    assert movies != []
