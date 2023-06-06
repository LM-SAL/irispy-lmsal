import pytest

from irispy.data.test import get_test_filepath


@pytest.fixture()
def raster_file():
    return get_test_filepath("iris_l2_20210905_001833_3620258102_raster_t000_r00000_test.fits")


@pytest.fixture()
def sji_1330_file():
    return get_test_filepath("iris_l2_20210905_001833_3620258102_SJI_1330_t000_test.fits")


@pytest.fixture()
def sji_1400_file():
    return get_test_filepath("iris_l2_20210905_001833_3620258102_SJI_1400_t000_test.fits")


@pytest.fixture()
def sji_2796_file():
    return get_test_filepath("iris_l2_20210905_001833_3620258102_SJI_2796_t000_test.fits")


@pytest.fixture()
def sji_2832_file():
    return get_test_filepath("iris_l2_20210905_001833_3620258102_SJI_2832_t000_test.fits")


@pytest.fixture()
def SJICube_1330():
    from irispy.io.sji import read_sji_lvl2

    return read_sji_lvl2(get_test_filepath("iris_l2_20210905_001833_3620258102_SJI_1330_t000_test.fits"))


@pytest.fixture()
def SJICube_1400():
    from irispy.io.sji import read_sji_lvl2

    return read_sji_lvl2(get_test_filepath("iris_l2_20210905_001833_3620258102_SJI_1400_t000_test.fits"))


@pytest.fixture()
def SJICube_2796():
    from irispy.io.sji import read_sji_lvl2

    return read_sji_lvl2(get_test_filepath("iris_l2_20210905_001833_3620258102_SJI_2796_t000_test.fits"))


@pytest.fixture()
def SJICube_2832():
    from irispy.io.sji import read_sji_lvl2

    return read_sji_lvl2(get_test_filepath("iris_l2_20210905_001833_3620258102_SJI_2832_t000_test.fits"))


@pytest.fixture()
def iris_responsev1():
    from sunpy.time import parse_time

    from irispy.utils import get_iris_response

    return get_iris_response(time_obs=parse_time("2013-09-03"), response_version=1)


@pytest.fixture()
def iris_responsev2():
    from sunpy.time import parse_time

    from irispy.utils import get_iris_response

    return get_iris_response(time_obs=parse_time("2013-09-03"), response_version=2)


@pytest.fixture()
def iris_responsev3():
    from sunpy.time import parse_time

    from irispy.utils import get_iris_response

    return get_iris_response(time_obs=parse_time("2013-09-03"), response_version=3)


@pytest.fixture()
def iris_responsev4():
    from sunpy.time import parse_time

    from irispy.utils import get_iris_response

    return get_iris_response(time_obs=parse_time("2013-09-03"), response_version=4)


@pytest.fixture()
def iris_responsev5():
    from sunpy.time import parse_time

    from irispy.utils import get_iris_response

    return get_iris_response(time_obs=parse_time("2013-09-03"), response_version=5)


@pytest.fixture()
def iris_responsev6():
    from sunpy.time import parse_time

    from irispy.utils import get_iris_response

    return get_iris_response(time_obs=parse_time("2013-09-03"), response_version=6)


@pytest.fixture()
def filelist():
    return [
        get_test_filepath("iris_l2_20210905_001833_3620258102_SJI_1330_t000_test.fits"),
        get_test_filepath("iris_l2_20210905_001833_3620258102_SJI_1400_t000_test.fits"),
        get_test_filepath("iris_l2_20210905_001833_3620258102_SJI_2796_t000_test.fits"),
        get_test_filepath("iris_l2_20210905_001833_3620258102_SJI_2832_t000_test.fits"),
    ]


@pytest.fixture(scope="session")
def fake_long_obs(tmp_path_factory):
    import os

    import numpy as np
    from astropy.io import fits

    header = fits.getheader(get_test_filepath("iris_l2_20210905_001833_3620258102_SJI_2832_t000_test.fits"))
    header["STARTOBS"] = "2017-05-02T05:25:51.000"
    header["ENDOBS"] = "2017-05-02T08:25:51.000"
    header["NAXIS3"] = 100
    rng = np.random.default_rng(12345)
    data = rng.random((header["NAXIS3"], header["NAXIS2"], header["NAXIS1"]))
    temp_dir = tmp_path_factory.mktemp("IRIS")
    hdu = fits.PrimaryHDU(data=data, header=header, do_not_scale_image_data=True, scale_back=True)
    fits_file = os.fspath(temp_dir.joinpath("iris_l2_20210905_001833_3620258102_SJI_2832_t000_test.fits"))
    hdu.writeto(fits_file)
    return [fits_file]
