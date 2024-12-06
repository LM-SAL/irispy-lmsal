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


def test_world_axis_physical_types_SJICube_2832(SJICube_2832):
    assert np.all(SJICube_2832.shape == (10, 40, 37))
    assert SJICube_2832.array_axis_physical_types == AXIS


def test_world_axis_physical_types_SJICube_2796(SJICube_2796):
    assert np.all(SJICube_2796.shape == (62, 40, 37))
    assert SJICube_2796.array_axis_physical_types == AXIS


def test_world_axis_physical_types_SJICube_1400(SJICube_1400):
    assert np.all(SJICube_1400.shape == (62, 40, 37))
    assert SJICube_1400.array_axis_physical_types == AXIS


def test_world_axis_physical_types_SJICube_1330(SJICube_1330):
    assert np.all(SJICube_1330.shape == (52, 40, 37))
    assert SJICube_1330.array_axis_physical_types == AXIS
