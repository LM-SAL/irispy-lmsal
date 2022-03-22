import pytest

from irispy.data.test import get_test_filepath
from irispy.io.utils import fitsinfo, read_files


@pytest.fixture
def raster_file():
    return get_test_filepath("iris_l2_20170502_052551_3893010094_raster_t000_r00000.fits")


@pytest.fixture
def sji_1330_file():
    return get_test_filepath("iris_l2_20170502_052551_3893010094_SJI_1330_t000.fits")


@pytest.fixture
def sji_1400_file():
    return get_test_filepath("iris_l2_20170502_052551_3893010094_SJI_1400_t000.fits")


@pytest.fixture
def sji_2796_file():
    return get_test_filepath("iris_l2_20170502_052551_3893010094_SJI_2796_t000.fits")


@pytest.fixture
def sji_2832_file():
    return get_test_filepath("iris_l2_20170502_052551_3893010094_SJI_2832_t000.fits")


def test_fitsinfo(capsys, raster_file, sji_1330_file, sji_1400_file, sji_2796_file, sji_2832_file):
    fitsinfo(raster_file)
    captured = capsys.readouterr()
    assert raster_file in captured.out

    fitsinfo(sji_1330_file)
    captured = capsys.readouterr()
    assert sji_1330_file in captured.out

    fitsinfo(sji_1400_file)
    captured = capsys.readouterr()
    assert sji_1400_file in captured.out

    fitsinfo(sji_2796_file)
    captured = capsys.readouterr()
    assert sji_2796_file in captured.out

    fitsinfo(sji_2832_file)
    captured = capsys.readouterr()
    assert sji_2832_file in captured.out


def test_read_files_errors_with_mix(raster_file, sji_1330_file):
    with pytest.raises(ValueError, match="You cannot mix raster and SJI files."):
        read_files([raster_file, sji_1330_file])


def test_read_files_raster(raster_file):
    # Simple test to ensure it doesnt error
    read_files(raster_file)
    read_files([raster_file])


def test_read_files_sji(sji_1330_file, sji_1400_file, sji_2796_file, sji_2832_file):
    # Simple test to ensure it doesnt error
    read_files(sji_1330_file)
    read_files(sji_1400_file)
    read_files(sji_2796_file)
    read_files(sji_2832_file)
    read_files([sji_2832_file])


def test_read_files_sji_more_than_one(sji_1330_file):
    # Don't allow more than one SJI file for now.
    with pytest.raises(ValueError, match="You cannot load more than one SJI file at a time."):
        read_files([sji_1330_file, sji_1330_file])
