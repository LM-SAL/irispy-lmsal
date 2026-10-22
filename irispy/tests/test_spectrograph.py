import io
import re
import copy

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pytest

import astropy.units as u
from astropy.io import fits
from astropy.tests.helper import assert_quantity_allclose
from astropy.wcs import WCS

from ndcube import NDCubeSequence

from irispy.io.spectrograph import read_spectrograph_lvl2
from irispy.io.utils import read_files
from irispy.spectrograph import SpectrogramCube, SpectrogramCubeSequence
from irispy.tests.helpers import figure_test, make_test_spectrogram_cube
from irispy.utils.constants import SLIT_WIDTH
from irispy.visualization import _registered_colormap


def test_fits_data_comparison(sns_sg_file):
    """
    Make sure the data is the same in pyfits and irispy.
    """
    iris_l2_test_raster = read_spectrograph_lvl2(sns_sg_file)
    with fits.open(sns_sg_file) as hdulist:
        spectral_window1 = hdulist[0].header["TDESC1"]
        spectral_window2 = hdulist[0].header["TDESC2"]
        spectral_window3 = hdulist[0].header["TDESC3"]
        data1 = copy.deepcopy(hdulist[1].data)
        data2 = copy.deepcopy(hdulist[2].data)
        data3 = copy.deepcopy(hdulist[3].data)
        np.testing.assert_array_almost_equal(iris_l2_test_raster[spectral_window1].data[0].data, data1)
        np.testing.assert_array_almost_equal(iris_l2_test_raster[spectral_window2].data[0].data, data2)
        np.testing.assert_array_almost_equal(iris_l2_test_raster[spectral_window3].data[0].data, data3)


def test_spectrogram_cube_slice_preserves_coordinates(sns_sg_file):
    raster = read_spectrograph_lvl2(sns_sg_file)
    cube = raster["C II 1336"][0]

    sliced_cube = cube[0]
    sliced_sequence = cube[:1]

    assert "time" in tuple(sliced_cube.global_coords)
    assert "time" in tuple(sliced_sequence.extra_coords.keys())


def test_time_is_the_time_coordinate(sns_sg_file):
    for (cube,) in read_files(sns_sg_file).values():
        for piece in (cube, cube[:, 5], cube[1:3, 4:9, 2:20], cube[0:1]):
            time, expected = piece.time, piece.axis_world_coords("time", wcs=piece.extra_coords)[0]
            assert (time.shape, time.scale, time.format) == (expected.shape, expected.scale, expected.format)
            # The same instants; the WCS splits them into jd1 and jd2 differently, so == can be off by a rounding
            assert_quantity_allclose((time - expected).to(u.s), 0 * u.s, atol=1 * u.ns)
        exposure, expected = cube[-1].time, cube.time[-1]
        assert (exposure.shape, exposure.scale, exposure.format) == ((), expected.scale, expected.format)
        assert exposure == expected


def test_spectrogram_cube_negative_indices(sns_sg_file):
    cube = read_spectrograph_lvl2(sns_sg_file)["C II 1336"][0]
    last = cube.shape[0] - 1

    def slit_bottom(sliced):
        return sliced.wcs.array_index_to_world(0, 0)[1]

    assert slit_bottom(cube[-1]) == slit_bottom(cube[last])
    assert slit_bottom(cube[1:5][-1]) == slit_bottom(cube[4])
    assert cube[-3:].shape == (3, *cube.shape[1:])
    with pytest.raises(IndexError, match="length-0 axis"):
        cube[3:3]
    for item in (slice(None, None, 2), slice(None, None, -1), (slice(None), slice(5, 1, -1))):
        with pytest.raises(IndexError, match="with a step is not supported"):
            cube[item]


def test_spectral_dispersion():
    wavelengths = np.arange(10) * 0.02 * u.nm + 140 * u.nm
    cube = make_test_spectrogram_cube(np.ones((1, 1, 10)), wavelengths)
    dispersion = cube.spectral_dispersion
    assert dispersion.unit.is_equivalent(u.nm)
    assert u.isclose(dispersion, 0.02 * u.nm, rtol=0.01)


def test_solid_angle():
    wavelengths = np.arange(10) * 0.02 * u.nm + 140 * u.nm
    cube = make_test_spectrogram_cube(np.ones((1, 1, 10)), wavelengths)
    angle = cube.solid_angle
    assert angle.unit.is_equivalent(u.sr)
    expected = 1.0 * u.arcsec * SLIT_WIDTH
    assert u.isclose(angle.to(u.sr), expected.to(u.sr), rtol=0.01)


def test_wavelength_axis():
    wavelengths = np.arange(10) * 0.02 * u.nm + 140 * u.nm
    cube = make_test_spectrogram_cube(np.ones((1, 1, 10)), wavelengths)
    assert cube.wavelength_axis == 2


def test_spectrogram_sequence_slices_as_ndcube_sequence(sns_sg_file):
    # Any slice but an integer raised the OBSID check on the empty sequence ndcube starts from
    cube = read_spectrograph_lvl2(sns_sg_file)["C II 1336"][0]
    sequence = SpectrogramCubeSequence([cube, cube], common_axis=0)
    reference = NDCubeSequence([cube, cube], common_axis=0)
    for sliced, expected in (
        (sequence[0:2], reference[0:2]),
        (sequence[:, 0:5], reference[:, 0:5]),
        (sequence.index_as_cube[3:9], reference.index_as_cube[3:9]),
    ):
        assert isinstance(sliced, SpectrogramCubeSequence)
        assert [c.shape for c in sliced.data] == [c.shape for c in expected.data]


def test_length_one_step_axis(sns_sg_file):
    # The combined WCS of a length-1 step axis raises in ndcube
    cube = read_spectrograph_lvl2(sns_sg_file)["C II 1336"][0][0:1]
    assert cube.wavelength_axis == 2
    assert "('em.wl',)" in str(cube)


@pytest.mark.parametrize(
    ("property_name", "message"),
    [("wavelength_axis", "wavelength axis"), ("spectral_dispersion", "no WAVE ctype")],
)
def test_spectral_properties_raise_without_wave(property_name, message):
    header = fits.Header()
    header["NAXIS"] = 2
    header["NAXIS1"] = 5
    header["NAXIS2"] = 2
    header["CTYPE1"] = "HPLT-TAN"
    header["CTYPE2"] = "HPLN-TAN"
    header["CDELT1"] = 0.1
    header["CRVAL1"] = 0
    header["CRPIX1"] = 1
    header["CUNIT1"] = "arcsec"
    header["CDELT2"] = 0.1
    header["CRVAL2"] = 0
    header["CRPIX2"] = 1
    header["CUNIT2"] = "arcsec"
    wcs = WCS(header)
    cube = SpectrogramCube(np.ones((2, 5)), wcs=wcs, uncertainty=None, unit=u.DN, meta={}, mask=None)
    with pytest.raises(ValueError, match=message):
        getattr(cube, property_name)


def test_solid_angle_missing_hplt_raises():
    header = fits.Header()
    header["NAXIS"] = 3
    header["NAXIS1"] = 5
    header["NAXIS2"] = 2
    header["NAXIS3"] = 1
    header["CTYPE1"] = "WAVE"
    header["CTYPE2"] = "TIME"
    header["CTYPE3"] = "UTC"
    header["CDELT1"] = 0.02
    header["CRVAL1"] = 140.0
    header["CRPIX1"] = 1
    header["CUNIT1"] = "nm"
    header["CDELT2"] = 1.0
    header["CRVAL2"] = 0
    header["CRPIX2"] = 1
    header["CUNIT2"] = "s"
    header["CDELT3"] = 1.0
    header["CRVAL3"] = 0
    header["CRPIX3"] = 1
    header["CUNIT3"] = "s"
    wcs = WCS(header)
    cube = SpectrogramCube(np.ones((1, 2, 5)), wcs=wcs, uncertainty=None, unit=u.DN, meta={}, mask=None)
    with pytest.raises(ValueError, match="no HPLT ctype"):
        _ = cube.solid_angle


def test_spectral_dispersion_real_data(sns_sg_file):
    raster = read_files(sns_sg_file)
    cube = raster["Si IV 1403"][0]
    dispersion = cube.spectral_dispersion
    assert dispersion.unit.is_equivalent(u.nm)
    assert dispersion.value > 0


def test_solid_angle_real_data(sns_sg_file):
    raster = read_files(sns_sg_file)
    cube = raster["Si IV 1403"][0]
    angle = cube.solid_angle
    assert angle.unit.is_equivalent(u.sr)
    assert angle.value > 0


def test_spectrogram_cube_plot_defaults(sns_sg_file):
    cube = read_files(sns_sg_file)["Si IV 1403"][0][0]
    # The IRIS FUV color table, as sunpy registers it, for an FUV window
    assert cube.meta.detector_band == "FUV"
    assert cube.plotter._default_cmap_name() == "irissjiFUV"
    ax = cube.plot()
    assert ax.images[0].get_cmap().name == "irissjiFUV"
    plt.close(ax.figure)


@pytest.mark.parametrize(
    ("kwargs", "labels"),
    [({}, ["283.20", "283.35"]), ({"aspect": "auto"}, ["283.20", "283.25", "283.30", "283.35", "283.40"])],
)
def test_spectrogram_plot_labels_narrow_window(sns_sg_file, kwargs, labels):
    # The 2832 window spans 0.23 nm, which the default equal aspect draws 92 pixels wide
    ax = read_files(sns_sg_file)["2832"][0][0].plot(**kwargs)
    svg = io.StringIO()
    with mpl.rc_context({"svg.fonttype": "none"}):
        ax.figure.savefig(svg, format="svg")
    plt.close(ax.figure)
    assert re.findall(r">(283\.\d+)</text>", svg.getvalue()) == labels


def test_default_colormap_falls_back_to_viridis():
    assert _registered_colormap("irissjiNUV") == "irissjiNUV"
    assert _registered_colormap("irissji1700") == "viridis"


@figure_test
def test_spectrogram_figure(sns_sg_file):
    # One exposure of the Si IV 1403 window: slit position against wavelength
    read_files(sns_sg_file)["Si IV 1403"][0][0].plot()


@pytest.mark.filterwarnings("ignore:Attempting to set identical low and high xlims")
def test_spectrogram_sequence_plot_uses_iris_slider_labels(sns_sg_file):
    sequence = read_files(sns_sg_file)["Si IV 1403"]
    fig = plt.figure()
    animator = sequence.plot(fig=fig)
    assert animator.slider_labels == ["Raster step", "Scan number"]
    plt.close(fig)


def test_raster_collection_aligned_axis_physical_types_are_sorted(sns_sg_file):
    # NDCollection builds these from sets, so without sorting the order changes between runs.
    types = read_files(sns_sg_file).aligned_axis_physical_types
    assert any(len(axis_types) > 1 for axis_types in types)
    assert types == [tuple(sorted(axis_types)) for axis_types in types]


def test_spectrogram_cube_from_data_and_wcs():
    # For example, a map of fitted line parameters with the WCS of the cube they came from.
    wcs = WCS(naxis=2)
    wcs.wcs.ctype = ["HPLT-TAN", "HPLN-TAN"]
    cube = SpectrogramCube(np.ones((4, 5)) * u.km / u.s, wcs)
    assert cube.unit == u.km / u.s
    assert cube.uncertainty is None
    assert str(cube)
