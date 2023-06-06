import copy

import astropy.units as u
import numpy as np
from astropy.io import fits
from astropy.time import Time, TimeDelta
from astropy.wcs import WCS

import irispy.data.test
from irispy import SpectrogramCube, SpectrogramCubeSequence, utils
from irispy.io.spectrograph import read_spectrograph_lvl2

testpath = irispy.data.test.rootdir
# Arrays of DN
SOURCE_DATA_DN = np.array(
    [
        [[0.563, 1.132, -1.343], [-0.719, 1.441, 1.566]],
        [[0.563, 1.132, -1.343], [-0.719, 1.441, 1.566]],
    ],
)
SOURCE_UNCERTAINTY_DN = np.sqrt(SOURCE_DATA_DN)
# Arrays relating SOURCE_DATA_DN to photons in NUV and FUV
SOURCE_DATA_PHOTONS_NUV = np.array(
    [
        [[10.134, 20.376, -24.174], [-12.942, 25.938, 28.188]],
        [[10.134, 20.376, -24.174], [-12.942, 25.938, 28.188]],
    ],
)
SOURCE_DATA_PHOTONS_FUV = np.array(
    [
        [[2.252, 4.528, -5.372], [-2.876, 5.764, 6.264]],
        [[2.252, 4.528, -5.372], [-2.876, 5.764, 6.264]],
    ],
)
SOURCE_UNCERTAINTY_PHOTONS_NUV = np.sqrt(SOURCE_DATA_PHOTONS_NUV)
SOURCE_UNCERTAINTY_PHOTONS_FUV = np.sqrt(SOURCE_DATA_PHOTONS_FUV)
time_dim_len = SOURCE_DATA_DN.shape[0]
single_exposure_time = 2.0
EXPOSURE_TIME = u.Quantity(np.zeros(time_dim_len) + single_exposure_time, unit=u.s)
# Define an sample wcs object
h0 = {
    "CTYPE1": "WAVE    ",
    "CUNIT1": "Angstrom",
    "CDELT1": 0.2,
    "CRPIX1": 0,
    "CRVAL1": 10,
    "NAXIS1": 3,
    "CTYPE2": "HPLT-TAN",
    "CUNIT2": "deg",
    "CDELT2": 0.5,
    "CRPIX2": 2,
    "CRVAL2": 0.5,
    "NAXIS2": 2,
    "CTYPE3": "HPLN-TAN",
    "CUNIT3": "deg",
    "CDELT3": 0.4,
    "CRPIX3": 2,
    "CRVAL3": 1,
    "NAXIS3": 2,
}
wcs0 = WCS(header=h0, naxis=3)
# Define sample meta
meta0 = {"detector type": "FUV", "OBSID": 1, "spectral window": "C II 1336"}
# Define sample extra coords
extra_coords0 = [
    ("time", 0, Time("2017-01-01") + TimeDelta(np.arange(time_dim_len), format="sec")),
    ("exposure time", 0, EXPOSURE_TIME),
]
extra_coords1 = [
    (
        "time",
        0,
        (Time("2017-01-01") + TimeDelta(np.arange(time_dim_len, time_dim_len * 2), format="sec")),
    ),
    ("exposure time", 0, EXPOSURE_TIME),
]
# Define SpectrogramCubes in various units.
spectrogram_DN0 = SpectrogramCube(
    SOURCE_DATA_DN,
    wcs0,
    SOURCE_UNCERTAINTY_DN,
    utils.DN_UNIT["FUV"],
    meta0,
)
spectrogram_DN0.extra_coords.add(*extra_coords0[0])
spectrogram_DN0.extra_coords.add(*extra_coords0[1])
spectrogram_photon0 = SpectrogramCube(
    SOURCE_DATA_PHOTONS_FUV,
    wcs0,
    SOURCE_UNCERTAINTY_PHOTONS_FUV,
    u.photon,
    meta0,
)
spectrogram_DN0.extra_coords.add(*extra_coords0[0])
spectrogram_DN0.extra_coords.add(*extra_coords0[1])
spectrogram_DN_per_s0 = SpectrogramCube(
    SOURCE_DATA_DN / single_exposure_time,
    wcs0,
    SOURCE_UNCERTAINTY_DN / single_exposure_time,
    utils.DN_UNIT["FUV"] / u.s,
    meta0,
)
spectrogram_DN_per_s0.extra_coords.add(*extra_coords0[0])
spectrogram_DN_per_s0.extra_coords.add(*extra_coords0[1])
spectrogram_photon_per_s0 = SpectrogramCube(
    SOURCE_DATA_PHOTONS_FUV / single_exposure_time,
    wcs0,
    SOURCE_UNCERTAINTY_PHOTONS_FUV / single_exposure_time,
    u.photon / u.s,
    meta0,
)
spectrogram_photon_per_s0.extra_coords.add(*extra_coords0[0])
spectrogram_photon_per_s0.extra_coords.add(*extra_coords0[1])
spectrogram_DN1 = SpectrogramCube(
    SOURCE_DATA_DN,
    wcs0,
    SOURCE_UNCERTAINTY_DN,
    utils.DN_UNIT["FUV"],
    meta0,
)
spectrogram_DN1.extra_coords.add(*extra_coords1[0])
spectrogram_DN1.extra_coords.add(*extra_coords1[1])
spectrogram_photon1 = SpectrogramCube(
    SOURCE_DATA_PHOTONS_FUV,
    wcs0,
    SOURCE_UNCERTAINTY_PHOTONS_FUV,
    u.photon,
    meta0,
)
spectrogram_photon1.extra_coords.add(*extra_coords1[0])
spectrogram_photon1.extra_coords.add(*extra_coords1[1])
spectrogram_DN_per_s1 = SpectrogramCube(
    SOURCE_DATA_DN / single_exposure_time,
    wcs0,
    SOURCE_UNCERTAINTY_DN / single_exposure_time,
    utils.DN_UNIT["FUV"] / u.s,
    meta0,
)
spectrogram_DN_per_s1.extra_coords.add(*extra_coords1[0])
spectrogram_DN_per_s1.extra_coords.add(*extra_coords1[1])
spectrogram_photon_per_s1 = SpectrogramCube(
    SOURCE_DATA_PHOTONS_FUV / single_exposure_time,
    wcs0,
    SOURCE_UNCERTAINTY_PHOTONS_FUV / single_exposure_time,
    u.photon / u.s,
    meta0,
)
spectrogram_DN0.extra_coords.add(*extra_coords0[0])
spectrogram_DN0.extra_coords.add(*extra_coords0[1])
spectrogram_photon_per_s_per_s0 = SpectrogramCube(
    SOURCE_DATA_PHOTONS_FUV / single_exposure_time / single_exposure_time,
    wcs0,
    SOURCE_UNCERTAINTY_PHOTONS_FUV / single_exposure_time / single_exposure_time,
    u.photon / u.s / u.s,
    meta0,
)
spectrogram_photon_per_s_per_s0.extra_coords.add(*extra_coords0[0])
spectrogram_photon_per_s_per_s0.extra_coords.add(*extra_coords0[1])
spectrogram_photon_s0 = SpectrogramCube(
    SOURCE_DATA_PHOTONS_FUV * single_exposure_time,
    wcs0,
    SOURCE_UNCERTAINTY_PHOTONS_FUV * single_exposure_time,
    u.photon * u.s,
    meta0,
)
spectrogram_photon_s0.extra_coords.add(*extra_coords0[0])
spectrogram_photon_s0.extra_coords.add(*extra_coords0[1])
spectrogram_photon_per_s_per_s1 = SpectrogramCube(
    SOURCE_DATA_PHOTONS_FUV / single_exposure_time / single_exposure_time,
    wcs0,
    SOURCE_UNCERTAINTY_PHOTONS_FUV / single_exposure_time / single_exposure_time,
    u.photon / u.s / u.s,
    meta0,
)
spectrogram_photon_per_s_per_s1.extra_coords.add(*extra_coords1[0])
spectrogram_photon_per_s_per_s1.extra_coords.add(*extra_coords1[1])
spectrogram_photon_s1 = SpectrogramCube(
    SOURCE_DATA_PHOTONS_FUV * single_exposure_time,
    wcs0,
    SOURCE_UNCERTAINTY_PHOTONS_FUV * single_exposure_time,
    u.photon * u.s,
    meta0,
)
spectrogram_photon_s1.extra_coords.add(*extra_coords1[0])
spectrogram_photon_s1.extra_coords.add(*extra_coords1[1])
# Define meta dict for an SpectrogramCubeSequence
meta_seq = {
    "detector type": "FUV",
    "spectral window": "C II 1336",
    "brightest wavelength": 100,
    "min wavelength": 90,
    "max wavelength": 110,
}
# Define SpectrogramCubeSequences
sequence_DN = SpectrogramCubeSequence([spectrogram_DN0, spectrogram_DN1], meta_seq)
sequence_photon = SpectrogramCubeSequence([spectrogram_photon0, spectrogram_photon1], meta_seq)
sequence_DN_per_s = SpectrogramCubeSequence([spectrogram_DN_per_s0, spectrogram_DN_per_s1], meta_seq)
sequence_photon_per_s = SpectrogramCubeSequence([spectrogram_photon_per_s0, spectrogram_photon_per_s1], meta_seq)
sequence_photon_per_s_per_s = SpectrogramCubeSequence(
    [spectrogram_photon_per_s_per_s0, spectrogram_photon_per_s1],
    meta_seq,
)
sequence_photon_s = SpectrogramCubeSequence([spectrogram_photon_s0, spectrogram_photon_s1], meta_seq)


def test_fits_data_comparison(raster_file):
    """
    Make sure the data is the same in pyfits and irispy.
    """
    iris_l2_test_raster = read_spectrograph_lvl2(raster_file)
    with fits.open(raster_file) as hdulist:
        spectral_window1 = hdulist[0].header["TDESC1"]
        spectral_window2 = hdulist[0].header["TDESC2"]
        spectral_window3 = hdulist[0].header["TDESC3"]
        data1 = copy.deepcopy(hdulist[1].data)
        data2 = copy.deepcopy(hdulist[2].data)
        data3 = copy.deepcopy(hdulist[3].data)
        np.testing.assert_array_almost_equal(iris_l2_test_raster[spectral_window1].data[0].data, data1)
        np.testing.assert_array_almost_equal(iris_l2_test_raster[spectral_window2].data[0].data, data2)
        np.testing.assert_array_almost_equal(iris_l2_test_raster[spectral_window3].data[0].data, data3)
