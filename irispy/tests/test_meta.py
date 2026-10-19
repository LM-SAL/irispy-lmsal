import numpy as np
import pytest

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.time import Time

from sunpy.coordinates import Helioprojective, get_earth
from sunpy.coordinates.utils import get_heliocentric_angle

from irispy.io.utils import read_files
from irispy.meta import SGMeta, SJIMeta


def _make_sg_header():
    header = fits.Header()
    header["INSTRUME"] = "IRIS"
    header["TELESCOP"] = "IRIS"
    header["DATA_LEV"] = 2
    header["OBSID"] = 1
    header["OBS_DESC"] = "Test obs"
    header["NWIN"] = 1
    header["TDESC1"] = "Si IV 1403"
    header["TWAVE1"] = 1402.77
    header["TWMIN1"] = 1398.61
    header["TWMAX1"] = 1406.03
    header["TDET1"] = "FUV2"
    return header


def test_sgmeta_rest_wavelength():
    meta = SGMeta(_make_sg_header(), "Si IV 1403")
    assert meta.rest_wavelength.unit == u.nm
    assert u.isclose(meta.rest_wavelength, 140.277 * u.nm, rtol=1e-4)


@pytest.mark.parametrize(
    ("spectral_window", "rest_wavelength"),
    [("Mg II k 2796", 2796.352 * u.AA), ("2832", 2832.982 * u.AA)],
    ids=["line", "continuum"],
)
def test_sgmeta_rest_wavelength_from_line_catalog(sns_sg_file, spectral_window, rest_wavelength):
    meta = SGMeta(fits.getheader(sns_sg_file), spectral_window)
    assert u.isclose(meta.rest_wavelength, rest_wavelength, atol=1 * u.mAA)


def test_sgmeta_detector_band_fuv():
    meta = SGMeta(_make_sg_header(), "Si IV 1403")
    assert meta.detector_band == "FUV"


def test_sgmeta_detector_band_nuv():
    header = _make_sg_header()
    header["TDET1"] = "NUV"
    meta = SGMeta(header, "Si IV 1403")
    assert meta.detector_band == "NUV"


def test_sgmeta_detector_band_sji():
    header = fits.Header()
    header["INSTRUME"] = "IRIS"
    header["TELESCOP"] = "IRIS"
    header["TDET1"] = "SJI"
    header["TDESC1"] = "SJI_1330"
    header["TWAVE1"] = 1330.0
    header["TWMIN1"] = 1330.0
    header["TWMAX1"] = 1330.0
    meta = SJIMeta(header)
    assert meta.detector_band == "SJI"


@pytest.mark.parametrize(
    ("spectral_window", "message"),
    [("Mg II k 2796", "Spectral window not found"), ("Si IV", "Spectral window must be unique")],
    ids=["missing", "ambiguous"],
)
def test_sgmeta_rejects_missing_or_ambiguous_window(spectral_window, message):
    header = _make_sg_header()
    header["NWIN"] = 2
    header["TDESC2"] = "Si IV 1394"
    with pytest.raises(ValueError, match=message):
        SGMeta(header, spectral_window)


def test_sgmeta_observer_uses_campaign_start_when_date_obs_is_blank():
    header = _make_sg_header()
    header["DATE_OBS"] = ""
    header["STARTOBS"] = "2014-03-29T14:09:38.830"
    header["XCEN"] = 0.0
    header["YCEN"] = 0.0
    meta = SGMeta(header, "Si IV 1403")

    assert meta.date_reference is None
    assert meta.observer.obstime.isot == "2014-03-29T14:09:38.830"


def test_sgmeta_real_sns_data(sns_sg_file):
    raster = read_files(sns_sg_file)
    cube = raster["Si IV 1403"][0]
    meta = cube.meta

    assert meta.detector_band == "FUV"
    assert meta.rest_wavelength.unit == u.nm
    assert u.isclose(meta.rest_wavelength, 140.277 * u.nm, rtol=1e-3)
    assert meta.temporal_cadence.unit.is_equivalent(u.s)
    assert meta.camera == 1
    assert meta.number_of_spectral_windows is not None


def test_sgmeta_sun_angular_radius_from_dsun():
    header = _make_sg_header()
    header["DSUN_OBS"] = 1.5e11
    meta = SGMeta(header, "Si IV 1403", data_shape=(10, 2, 2))
    radius = meta.sun_angular_radius
    assert radius.unit.is_equivalent(u.arcsec)


@pytest.mark.parametrize(
    ("xcen", "ycen", "expected"),
    [(0, 0, 1), (0, 960, 0), (-1000, 0, np.nan)],
    ids=["disk center", "limb", "off disk"],
)
def test_meta_mu(xcen, ycen, expected):
    header = _make_sg_header()
    header["RSUN_OBS"] = 960.0
    header["XCEN"], header["YCEN"] = xcen, ycen
    np.testing.assert_allclose(SGMeta(header, "Si IV 1403").mu, expected, atol=1e-12)


@pytest.mark.parametrize("fraction", [0.5, 0.999])
def test_meta_mu_matches_sunpy(fraction):
    observer = get_earth("2021-10-01T06:09:25")
    header = _make_sg_header()
    header["DSUN_OBS"] = observer.radius.to_value(u.m)
    # On the diagonal, where the angular distance from disk center is not sqrt(Tx**2 + Ty**2)
    offset = fraction * SGMeta(header, "Si IV 1403").sun_angular_radius / np.sqrt(2)
    header["XCEN"] = header["YCEN"] = offset.to_value(u.arcsec)
    point = SkyCoord(offset, offset, frame=Helioprojective(observer=observer, obstime=observer.obstime))
    expected = np.cos(get_heliocentric_angle(point)).to_value(u.one)
    np.testing.assert_allclose(SGMeta(header, "Si IV 1403").mu, expected, atol=1e-8)


def test_sgmeta_exposure_mu(raster_sg_file):
    cube = read_files(raster_sg_file)["C II 1336"][0]
    exposure_mu = cube.meta.exposure_mu
    assert exposure_mu.shape == cube.data.shape[:1]
    # The raster steps across its center, so mu at the center lies within the per-step range
    assert exposure_mu.min() < cube.meta.mu < exposure_mu.max()
    np.testing.assert_array_equal(cube[1:3].meta.exposure_mu, exposure_mu[1:3])


def test_sgmeta_observer_radial_velocity_per_exposure(sns_sg_file):
    # Level 2 primary headers have no OBS_VR, only the per-exposure OBS_VRIX
    meta = read_files(sns_sg_file)["Si IV 1403"][0].meta
    with fits.open(sns_sg_file) as hdulist:
        assert "OBS_VR" not in hdulist[0].header
        expected = hdulist[-2].data[:, hdulist[-2].header["OBS_VRIX"]] * u.m / u.s
    assert u.allclose(meta.observer_radial_velocity, expected)


@pytest.mark.parametrize(
    ("keyword", "value", "name", "expected"),
    [
        ("CAMERA", 1, "camera", 1),
        ("RSUN_OBS", 975.0, "sun_angular_radius", 975.0 * u.arcsec),
        ("CADEX_AV", 9.264, "temporal_cadence", 9.264 * u.s),
        ("OBS_VR", 3500.0, "observer_radial_velocity", 3500.0 * u.m / u.s),
        ("NWIN", 1, "number_of_spectral_windows", 1),
        ("RASRPT", 3, "raster_repetition", 3),
        ("EXPTIME", 12.5, "exposure_time", 12.5 * u.s),
        ("EXPMIN", 10.0, "exposure_time_min", 10.0 * u.s),
        ("EXPMAX", 15.0, "exposure_time_max", 15.0 * u.s),
        ("BTYPE", "Intensity", "data_type", "Intensity"),
        ("BUNIT", "DN/s", "data_unit", "DN/s"),
        ("STATUS", "Final", "data_status", "Final"),
        ("BLD_VERS", "v9.4", "build_version", "v9.4"),
        ("VER_RF2", "1.2.3", "reformat_version", "1.2.3"),
        ("DATE_RF2", "2024-01-15T10:30:00", "reformat_date", Time("2024-01-15T10:30:00", scale="utc")),
        ("OBSLABEL", "OBS_12345", "observing_label", "OBS_12345"),
        ("OBSTITLE", "Test Campaign", "observing_title", "Test Campaign"),
        ("LUTID", 7, "lut_id", 7),
        ("NEXP", 180, "number_of_exposures", 180),
        ("NEXP_PRP", 200, "number_of_exposures_planned", 200),
        ("NEXPOBS", 400, "number_of_exposures_observation", 400),
        ("NSATPIX", 42, "number_of_saturated_pixels", 42),
        ("NSPIKES", 15, "number_of_spikes", 15),
        ("PERCENTD", 98.5, "percent_data", 98.5),
        ("DATAMEAN", 120.5, "data_mean", 120.5),
        ("DATARMS", 45.2, "data_rms", 45.2),
        ("DATAMEDN", 110.0, "data_median", 110.0),
        ("DATAMIN", 5.0, "data_min", 5.0),
        ("DATAMAX", 500.0, "data_max", 500.0),
        ("STEPS_AV", 0.35, "step_size_average", 0.35),
        ("STEPS_DV", 0.02, "step_size_stddev", 0.02),
        ("STEPT_AV", 18.5, "step_time_average", 18.5),
        ("STEPT_DV", 1.2, "step_time_stddev", 1.2),
        ("CADPL_AV", 15.0, "cadence_planned_average", 15.0 * u.s),
        ("CADPL_DV", 0.5, "cadence_planned_stddev", 0.5 * u.s),
        ("CADEX_DV", 0.3, "cadence_executed_stddev", 0.3 * u.s),
        ("RASTYPDX", 3, "raster_type_index", 3),
        ("RASTYPNX", 8, "raster_type_total", 8),
        ("MISSRAS", 2, "number_of_missing_raster_files", 2),
        ("MISSOBS", 5, "number_of_missing_observation_files", 5),
        ("TDMEAN1", 150.0, "window_mean", 150.0),
        ("TDRMS1", 30.0, "window_rms", 30.0),
        ("TDMEDN1", 145.0, "window_median", 145.0),
        ("TDMIN1", 10.0, "window_min", 10.0),
        ("TDMAX1", 400.0, "window_max", 400.0),
        ("TSATPX1", 3, "window_saturated_pixels", 3),
        ("TSPIKE1", 7, "window_spikes", 7),
    ],
)
def test_sgmeta_reads_one_keyword(keyword, value, name, expected):
    header = _make_sg_header()
    header[keyword] = value
    meta = SGMeta(header, "Si IV 1403", data_shape=(10, 2, 2))

    result = getattr(meta, name)

    if isinstance(expected, u.Quantity):
        assert u.isclose(result, expected)
    elif isinstance(expected, Time):
        assert (result.isot, result.scale) == (expected.isot, expected.scale)
    else:
        assert result == expected


def test_sgmeta_missing_cadence_keys():
    meta = SGMeta(_make_sg_header(), "Si IV 1403", data_shape=(10, 2, 2))
    assert meta.temporal_cadence is None
    assert meta.cadence_planned_average is None
    assert meta.cadence_planned_stddev is None
    assert meta.cadence_executed_stddev is None
    assert meta.observer_radial_velocity is None


def test_sjimeta_aia_cutout_header():
    # The IRIS-aligned AIA cutouts lack TWMIN/TWMAX, CADEX_AV and OBS_VR, and store DATE_TIME_OBSID
    header = fits.Header()
    header["INSTRUME"] = "AIA_3"
    header["TDET1"] = "SJI"
    header["TDESC1"] = "171_THIN"
    header["TWAVE1"] = 171
    header["OBSID"] = "20250519_165924_3640107442"
    meta = SJIMeta(header, data_shape=(10, 10))
    assert meta.observing_mode_id == 3640107442
    for name in (
        "spectral_range",
        "temporal_cadence",
        "observer_radial_velocity",
        "distance_to_sun",
        "sun_angular_radius",
        "mu",
        "satellite_rotation",
        "exposure_time",
        "processing_level",
        "raster_fov_width_x",
        "raster_fov_width_y",
    ):
        assert getattr(meta, name) is None, name
    assert "Spectral Range:  None" in str(meta)


def test_sjmeta_base_properties():
    header = fits.Header()
    header["INSTRUME"] = "IRIS"
    header["TELESCOP"] = "IRIS"
    header["DATA_LEV"] = 2
    header["TDET1"] = "SJI"
    header["TDESC1"] = "SJI_1330"
    header["TWAVE1"] = 1330.0
    header["TWMIN1"] = 1330.0
    header["TWMAX1"] = 1330.0
    header["EXPTIME"] = 8.0
    header["BTYPE"] = "Intensity"
    header["NSATPIX"] = 0
    meta = SJIMeta(header, data_shape=(10, 10))
    assert u.isclose(meta.exposure_time, 8.0 * u.s)
    assert meta.data_type == "Intensity"
    assert meta.number_of_saturated_pixels == 0
