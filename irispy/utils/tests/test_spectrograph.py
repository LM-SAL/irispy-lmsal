import numpy as np
import numpy.testing as np_test
import pytest

import astropy.units as u
from astropy.tests.helper import assert_quantity_allclose
from irispy.io.utils import read_files
from irispy.utils.constants import DN_UNIT, SLIT_WIDTH
from irispy.utils.response import get_latest_response
from irispy.utils.spectrograph import calculate_dn_to_radiance_factor
from irispy.data.test import get_test_filepath
from scipy.io import readsav
from irispy.utils.utils import record_to_dict
from irispy.utils.constants import DN_UNIT


# def test_convert_from_dn_to_photons(raster_file):
#     raster_collection = read_files(raster_file)
#     cube_sequence = raster_collection["C II 1336"]
#     cube = cube_sequence[0]
#     raster_collection_ph = convert_from_dn_to_photons(raster_collection)
#     cube_sequence_ph = convert_from_dn_to_photons(cube_sequence)
#     cube_ph = convert_from_dn_to_photons(cube)

#     assert raster_collection_ph["C II 1336"].unit == u.photon
#     assert cube_sequence_ph[0].unit == u.photon
#     assert cube_ph.unit == u.photon

# def test_radiometric_calibration(raster_file):
#     raster_cube = read_files(raster_file)
#     # Fix this later.
#     cali_results = radiometric_calibration(raster_cube["C II 1336"][0])
#     assert cali_results is not None


@pytest.fixture
def idl_input_rad_cal():
    return readsav(get_test_filepath("input_calibration.sav"))


@pytest.fixture
def idl_output_rad_cal():
    return readsav(get_test_filepath("output_calibration.sav"))

@pytest.fixture
def idl_output_factor_si():
    return readsav(get_test_filepath("output_factor.sav"))


@pytest.fixture
def idl_flux_radcal():
    return readsav(get_test_filepath("iris_datacal.sav"))


# def test_convert_photons_per_sec_to_radiance(raster_file, idl_input_rad_cal, idl_output_rad_cal):
#     raster_collection = read_files(raster_file)
#     cube = raster_collection["C II 1336"][0]
#     breakpoint()
#     idl_input_spectrum, idl_wavelength = idl_input_rad_cal["input_spectrum"], idl_input_rad_cal["wavelength"] * u.Angstrom
#     idl_calib_spectrum = idl_output_rad_cal["calib_spectrum"]

#     spectral_wcs_index = np.where(np.array(cube.wcs.wcs.ctype) == "WAVE")[0][0]
#     spectral_dispersion_per_pixel = cube.wcs.wcs.cdelt[spectral_wcs_index] * cube.wcs.wcs.cunit[spectral_wcs_index]
#     # Get solid angle from slit width for a pixel.
#     lat_wcs_index = ["HPLT" in c for c in cube.wcs.wcs.ctype]
#     lat_wcs_index = np.arange(len(cube.wcs.wcs.ctype))[lat_wcs_index]
#     lat_wcs_index = lat_wcs_index[0]
#     solid_angle = cube.wcs.wcs.cdelt[lat_wcs_index] * cube.wcs.wcs.cunit[lat_wcs_index] * SLIT_WIDTH
#     # Get wavelength for each pixel.
#     obs_wavelength = cube.axis_world_coords(2)
#     time_obs = cube.meta.date_reference

#     # Call the function to test
#     output = convert_photons_per_sec_to_radiance(data_quantities=(idl_input_spectrum * u.photon / u.s,),time_obs=time_obs,obs_wavelength=idl_wavelength,detector_type="FUV",spectral_dispersion_per_pixel=spectral_dispersion_per_pixel,solid_angle=solid_angle)

#     # Compare the output with the expected output
#     np_test.assert_allclose(output, expected_output)



#@pytest.mark.parametrize("window", ("A", "B"))
def test_calculate_dn_to_radiance_factor(raster_file, idl_input_rad_cal, idl_output_factor_si):
    raster_collection = read_files(raster_file)
    cube = raster_collection["C II 1336"][0]
    idl_wavelength = idl_input_rad_cal["wavelength"] * u.Angstrom
    idl_factor_si = idl_output_factor_si["factor"]

    spectral_wcs_index = np.where(np.array(cube.wcs.wcs.ctype) == "WAVE")[0][0]
    spectral_dispersion_per_pixel = cube.wcs.wcs.cdelt[spectral_wcs_index] * cube.wcs.wcs.cunit[spectral_wcs_index]

    # Get solid angle from slit width for a pixel.
    lat_wcs_index = ["HPLT" in c for c in cube.wcs.wcs.ctype]
    lat_wcs_index = np.arange(len(cube.wcs.wcs.ctype))[lat_wcs_index]
    lat_wcs_index = lat_wcs_index[0]
    solid_angle = cube.wcs.wcs.cdelt[lat_wcs_index] * cube.wcs.wcs.cunit[lat_wcs_index] * SLIT_WIDTH
    # These are using the date obs of now and not of the obs, i forgot to
    # pass that into the IDL function
    iris_response = get_latest_response(None)
    factor = calculate_dn_to_radiance_factor(
        iris_response,
        idl_wavelength,
        "FUV",
        spectral_dispersion_per_pixel,
        solid_angle,
    )
    assert len(factor) == len(idl_wavelength)
    # Python factor assumes data is in photons / s
    factor = factor * 4 * (u.photon / u.s)
    factor = factor.to("W / (m2 sr nm)")
    # IDL returns it in full units
    # W m^-2 s^-1 sr^-1 nm^-1.
    idl_factor = idl_factor_si * (u.W / (u.m**2 * u.sr * u.nm))
    # IRIS_CALIB_SPECTRUM: Spectrum converted by: (min,mean,max)=(28.811732,30.982474,33.369611)
    assert_quantity_allclose(factor, idl_factor)