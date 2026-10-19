import matplotlib.pyplot as plt
import numpy as np
import pytest

from astropy import units as u
from astropy.time import Time
from astropy.wcs import WCS

from ndcube.meta import NDMeta

from irispy import SJICube, utils
from irispy.utils.constants import BAD_PIXEL_VALUE_SCALED

TIMES = Time(["2014-12-11T19:39:00.48", "2014-12-11T19:43:07.6"])


def _fake_sji_cube(data, *, uncertainty, mask, **meta):
    """
    A slit-jaw cube on a plain FITS WCS, with the exposure times and time of two frames.
    """
    data = np.asarray(data)
    header = {
        "CTYPE1": "HPLN-TAN",
        "CUNIT1": "arcsec",
        "CDELT1": 0.4,
        "CRPIX1": 0,
        "CRVAL1": 0,
        "NAXIS1": data.shape[-1],
        "CTYPE2": "HPLT-TAN",
        "CUNIT2": "arcsec",
        "CDELT2": 0.5,
        "CRPIX2": 0,
        "CRVAL2": 0,
        "NAXIS2": data.shape[-2],
    }
    if data.ndim == 3:
        header |= {
            "CTYPE3": "Time    ",
            "CUNIT3": "s",
            "CDELT3": 0.3,
            "CRPIX3": 0,
            "CRVAL3": 0,
            "NAXIS3": data.shape[0],
        }
    cube = SJICube(
        data,
        WCS(header=header, naxis=data.ndim),
        uncertainty=uncertainty,
        mask=mask,
        unit=utils.constants.DN_UNIT["SJI"],
        meta=NDMeta(
            {"exposure time": 2 * np.ones(2) * u.s, "scaled": True, **meta},
            axes={"exposure time": 0},
            data_shape=data.shape,
        ),
    )
    cube.extra_coords.add("time", 0, TIMES)
    return cube


@pytest.fixture
def cube():
    data = np.array(
        [
            [[1, 2, 3, 4], [2, 4, 5, 3], [0, 1, 2, 3]],
            [[2, 4, 5, 1], [10, 5, 2, 2], [10, 3, 3, 0]],
        ],
    )
    return _fake_sji_cube(data, uncertainty=np.sqrt(data), mask=data >= 0)


@pytest.fixture
def cube_2d():
    data = np.array([[1, 2, 3, 4], [2, 4, 5, 3]])
    return _fake_sji_cube(data, uncertainty=np.sqrt(data), mask=data >= 0)


@pytest.fixture
def dust_cube():
    data = np.array(
        [
            [[-1, 2, -3, 4], [2, BAD_PIXEL_VALUE_SCALED, 5, 3], [0, 1, 2, -300]],
            [[2, BAD_PIXEL_VALUE_SCALED, 5, 1], [10, -5, 2, 2], [10, -3, 3, 0]],
        ],
    )
    return _fake_sji_cube(data, uncertainty=1, mask=data == BAD_PIXEL_VALUE_SCALED, OBSID=1)


def test_sjicube_apply_dust_mask(dust_cube):
    # TODO: The expected values are not correct.
    dust_mask_expected = np.array(
        [
            [[True, True, True, True], [True, True, True, True], [True, True, False, False]],
            [[True, True, True, False], [True, True, True, True], [True, True, True, True]],
        ]
    )
    dust_cube.apply_dust_mask()
    np.testing.assert_array_equal(dust_cube.mask, dust_mask_expected)
    dust_cube.apply_dust_mask(undo=True)
    before_mask = np.array(
        [
            [[False, False, False, False], [False, False, False, False], [False, False, False, False]],
            [[False, False, False, False], [False, False, False, False], [False, False, False, False]],
        ]
    )
    np.testing.assert_array_equal(dust_cube.mask, before_mask)


def test_sjicube_apply_dust_mask_initializes_missing_mask(dust_cube):
    dust_cube.mask = None

    dust_cube.apply_dust_mask()

    assert dust_cube.mask is not None
    assert dust_cube.mask.any()


@pytest.mark.parametrize(
    ("item", "expected_len"),
    [
        (-1, None),
        (0, None),
        (slice(0, 3), 3),
        ((slice(0, 3), slice(0, 10)), 3),
        ((slice(0, 3), slice(0, 10), slice(0, 10)), 3),
        ((0, slice(0, 10), slice(0, 10)), None),
        (Ellipsis, 52),
        ((Ellipsis, slice(0, 10)), 52),
        ((slice(0, 3), Ellipsis), 3),
        ((0, Ellipsis), None),
    ],
)
def test_sjicube_slice_preserves_fits_wcs(sns_sjicube_1330, item, expected_len):
    fits_wcs = sns_sjicube_1330.fits_wcs
    assert sns_sjicube_1330.fits_wcs is fits_wcs

    subset = sns_sjicube_1330[item]

    subset_fits_wcs = subset.fits_wcs
    assert subset_fits_wcs is not None
    assert subset.fits_wcs is subset_fits_wcs
    if expected_len is None:
        assert isinstance(subset.fits_wcs, WCS)
    else:
        assert len(subset.fits_wcs) == expected_len
        assert isinstance(subset.fits_wcs[0], WCS)


def test_sjicube_slice_rejects_multiple_ellipsis(sns_sjicube_1330):
    with pytest.raises((IndexError, ValueError), match=r"single ellipsis|only have a single ellipsis"):
        sns_sjicube_1330[(Ellipsis, Ellipsis)]


def test_sjicube_slice_rejects_step_slices(sns_sjicube_1330):
    with pytest.raises(IndexError, match=r"Slicing WCS with a step is not supported\."):
        sns_sjicube_1330[slice(0, 10, 2)]


def test_sjicube_slice_rejects_too_many_indices(sns_sjicube_1330):
    with pytest.raises((IndexError, ValueError), match=r"can not be greater than the dimensionality"):
        sns_sjicube_1330[(slice(0, 3), slice(0, 10), slice(0, 10), slice(None))]


def test_sjicube_slice_rejects_out_of_range_index(sns_sjicube_1330):
    with pytest.raises(IndexError, match=r"out of range for a cube of shape \(52, 40, 37\)") as error:
        sns_sjicube_1330[0, 100]
    assert not isinstance(error.value, ValueError)


def test_sjicube_slice_without_frame_headers_keeps_fits_wcs_none(cube):
    subset = cube[0:1]

    assert cube.fits_wcs is None
    assert subset.fits_wcs is None


def test_sjicube_plot_without_twave1_uses_viridis(cube, cube_2d):
    assert "TWAVE1" not in cube.meta  # fixture check
    assert cube.plotter._default_cmap_name() == "viridis"
    ax = cube_2d.plot()
    assert ax.images[0].get_cmap().name == "viridis"
    plt.close(ax.figure)


def test_sjicube_2d_slice(cube_2d):
    subset = cube_2d[:, :2]

    assert subset.data.ndim == 2
    assert subset.fits_wcs is None


def test_sjicube_remove_cosmic_rays(cube_2d, monkeypatch):
    captured = {}

    def fake_remove_cosmic_rays(data, mask, *, sigma, max_iters, method_kwargs):
        captured["mask"] = mask.copy()
        captured["sigma"] = sigma
        captured["max_iters"] = max_iters
        captured["method_kwargs"] = method_kwargs
        return data + 1

    monkeypatch.setattr("irispy.utils.cosmic_rays._remove_cosmic_rays_astroscrappy", fake_remove_cosmic_rays)

    original_data = cube_2d.data.copy()
    cleaned_cube = cube_2d.remove_cosmic_rays(
        method="astroscrappy",
        sigma=2.0,
        max_iters=3,
        method_kwargs={"readnoise": 4.0},
    )

    assert isinstance(cleaned_cube, SJICube)
    np.testing.assert_array_equal(cleaned_cube.data, cube_2d.data + 1)
    np.testing.assert_array_equal(cleaned_cube.mask, cube_2d.mask)
    assert cleaned_cube.meta["scaled"] == cube_2d.meta["scaled"]
    assert cleaned_cube.dust_masked == cube_2d.dust_masked
    assert list(cleaned_cube.extra_coords.keys()) == list(cube_2d.extra_coords.keys())
    assert captured["sigma"] == 2.0
    assert captured["max_iters"] == 3
    assert captured["method_kwargs"] == {"readnoise": 4.0}
    np.testing.assert_array_equal(captured["mask"], cube_2d.mask)
    np.testing.assert_array_equal(cube_2d.data, original_data)
    # WCS is deep-copied; check equivalence not identity
    assert dict(cleaned_cube.wcs.to_header()) == dict(cube_2d.wcs.to_header())
    assert cleaned_cube.fits_wcs == cube_2d.fits_wcs
    # Input mask must not be mutated
    np.testing.assert_array_equal(cube_2d.mask, original_data >= 0)
    # unit and meta keys are preserved on the cleaned cube
    assert cleaned_cube.unit == cube_2d.unit
    assert set(cleaned_cube.meta.keys()) == set(cube_2d.meta.keys())
