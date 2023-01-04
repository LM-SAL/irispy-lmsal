import astropy.units as u
import numpy as np

AXIS = [
    ("meta.obs.sequence",),
    ("custom:pos.helioprojective.lon", "custom:pos.helioprojective.lat"),
    ("custom:pos.helioprojective.lon", "custom:pos.helioprojective.lat"),
]


def test_world_axis_physical_types_IRISMapCube_2832(IRISMapCube_2832):
    assert np.all(IRISMapCube_2832.dimensions == [10, 40, 37] * u.pix)
    assert IRISMapCube_2832.array_axis_physical_types == AXIS


def test_world_axis_physical_types_IRISMapCube_2796(IRISMapCube_2796):
    assert np.all(IRISMapCube_2796.dimensions == [62, 40, 37] * u.pix)
    assert IRISMapCube_2796.array_axis_physical_types == AXIS


def test_world_axis_physical_types_IRISMapCube_1400(IRISMapCube_1400):
    assert np.all(IRISMapCube_1400.dimensions == [62, 40, 37] * u.pix)
    assert IRISMapCube_1400.array_axis_physical_types == AXIS


def test_world_axis_physical_types_IRISMapCube_1330(IRISMapCube_1330):
    assert np.all(IRISMapCube_1330.dimensions == [52, 40, 37] * u.pix)
    assert IRISMapCube_1330.array_axis_physical_types == AXIS
