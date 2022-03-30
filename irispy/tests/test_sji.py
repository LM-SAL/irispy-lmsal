import astropy.units as u
import numpy as np
import pytest

from irispy.data.test import get_test_filepath
from irispy.io.sji import read_sji_lvl2


@pytest.fixture
def IRISMapCube_1330():
    return read_sji_lvl2(get_test_filepath("iris_l2_20170502_052551_3893010094_SJI_1330_t000.fits"))


@pytest.fixture
def IRISMapCube_1400():
    return read_sji_lvl2(get_test_filepath("iris_l2_20170502_052551_3893010094_SJI_1400_t000.fits"))


@pytest.fixture
def IRISMapCube_2796():
    return read_sji_lvl2(get_test_filepath("iris_l2_20170502_052551_3893010094_SJI_2796_t000.fits"))


@pytest.fixture
def IRISMapCube_2832():
    return read_sji_lvl2(get_test_filepath("iris_l2_20170502_052551_3893010094_SJI_2832_t000.fits"))


AXIS = [
    ("meta.obs.sequence",),
    ("custom:pos.helioprojective.lon", "custom:pos.helioprojective.lat"),
    ("custom:pos.helioprojective.lon", "custom:pos.helioprojective.lat"),
]


def test_world_axis_physical_types_IRISMapCube_2832(IRISMapCube_2832):
    assert np.all(IRISMapCube_2832.dimensions == [16, 77, 148] * u.pix)
    assert IRISMapCube_2832.array_axis_physical_types == AXIS


def test_world_axis_physical_types_IRISMapCube_2796(IRISMapCube_2796):
    assert np.all(IRISMapCube_2796.dimensions == [16, 77, 148] * u.pix)
    assert IRISMapCube_2796.array_axis_physical_types == AXIS


def test_world_axis_physical_types_IRISMapCube_1400(IRISMapCube_1400):
    assert np.all(IRISMapCube_1400.dimensions == [16, 77, 148] * u.pix)
    assert IRISMapCube_1400.array_axis_physical_types == AXIS


def test_world_axis_physical_types_IRISMapCube_1330(IRISMapCube_1330):
    assert np.all(IRISMapCube_1330.dimensions == [16, 77, 148] * u.pix)
    assert IRISMapCube_1330.array_axis_physical_types == AXIS
