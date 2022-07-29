import pytest

from irispy.data.test import get_test_filepath


@pytest.fixture
def raster_file():
    return get_test_filepath("iris_l2_20210905_001833_3620258102_raster_t000_r00000.fits")


@pytest.fixture
def sji_1330_file():
    return get_test_filepath("iris_l2_20210905_001833_3620258102_SJI_1330_t000.fits")


@pytest.fixture
def sji_1400_file():
    return get_test_filepath("iris_l2_20210905_001833_3620258102_SJI_1400_t000.fits")


@pytest.fixture
def sji_2796_file():
    return get_test_filepath("iris_l2_20210905_001833_3620258102_SJI_2796_t000.fits")


@pytest.fixture
def sji_2832_file():
    return get_test_filepath("iris_l2_20210905_001833_3620258102_SJI_2832_t000.fits")
