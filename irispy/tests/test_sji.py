import numpy as np

AXIS = [
    (
        "custom:pos.helioprojective.lon",
        "custom:pos.helioprojective.lat",
        "time",
        "custom:CUSTOM",
        "custom:CUSTOM",
        "custom:CUSTOM",
        "custom:CUSTOM",
        "custom:CUSTOM",
        "custom:CUSTOM",
        "custom:CUSTOM",
        "custom:CUSTOM",
        "custom:CUSTOM",
    ),
    ("custom:pos.helioprojective.lon", "custom:pos.helioprojective.lat"),
    ("custom:pos.helioprojective.lon", "custom:pos.helioprojective.lat"),
]


def test_world_axis_physical_types_sjicube_2832(sjicube_2832):
    assert np.all(sjicube_2832.shape == (10, 40, 37))
    assert sjicube_2832.array_axis_physical_types == AXIS


def test_world_axis_physical_types_sjicube_2796(sjicube_2796):
    assert np.all(sjicube_2796.shape == (62, 40, 37))
    assert sjicube_2796.array_axis_physical_types == AXIS


def test_world_axis_physical_types_sjicube_1400(sjicube_1400):
    assert np.all(sjicube_1400.shape == (62, 40, 37))
    assert sjicube_1400.array_axis_physical_types == AXIS


def test_world_axis_physical_types_sjicube_1330(sjicube_1330):
    assert np.all(sjicube_1330.shape == (52, 40, 37))
    assert sjicube_1330.array_axis_physical_types == AXIS
