import os.path

import numpy as np
import numpy.testing as np_test
import pytest
import scipy.io
from sunpy.time import parse_time

from irispy import utils
from irispy.data.test import rootdir

# Arrays for the fit_iris_xput method
r = os.path.join(rootdir, "..", "iris_sra_c_20161022.geny")
raw_response_data = scipy.io.readsav(r, python_dict=True, verbose=False)
iris_response = dict([(name, raw_response_data["p0"][name][0]) for name in raw_response_data["p0"].dtype.names])
time_obs = parse_time("2013-09-03", format="utime")
time_cal_coeffs0 = iris_response.get("C_F_TIME")
time_cal_coeffs1 = iris_response.get("C_N_TIME")
cal_coeffs0 = iris_response.get("COEFFS_FUV")[0, :, :]
cal_coeffs1 = iris_response.get("COEFFS_FUV")[1, :, :]
cal_coeffs2 = iris_response.get("COEFFS_FUV")[2, :, :]
cal_coeffs3 = iris_response.get("COEFFS_NUV")[0, :, :]
cal_coeffs4 = iris_response.get("COEFFS_NUV")[1, :, :]
cal_coeffs5 = iris_response.get("COEFFS_NUV")[2, :, :]
iris_fit_expected0 = np.array([0.79865301])
iris_fit_expected1 = np.array([2.2495413])
iris_fit_expected2 = np.array([2.2495413])
iris_fit_expected3 = np.array([0.23529011])
iris_fit_expected4 = np.array([0.25203046])
iris_fit_expected5 = np.array([0.25265095])


def test_get_iris_response_not_equal_to_one():
    assert pytest.raises(
        KeyError,
        utils.get_iris_response,
        time_obs,
        pre_launch=False,
        response_version=13,
    )


def test_get_iris_response_response_file():
    assert pytest.raises(KeyError, utils.get_iris_response, time_obs, response_file="hello.py")


# Tests for get_iris_response function
# Version 1
sav_file_path1 = os.path.join(rootdir, "idl_iris_get_response_20130903_new_version001.sav")
test_iris_response1 = scipy.io.readsav(sav_file_path1, python_dict=True, verbose=False)
iris_response_load1 = test_iris_response1["iris_response"][0]

lamb_load1 = iris_response_load1.lambda_vars
area_sg_load1 = iris_response_load1.area_sg
name_sg_load1 = iris_response_load1.name_sg
index_el_sg_load1 = iris_response_load1.index_el_sg
area_sji_load1 = iris_response_load1.area_sji
name_sji_load1 = iris_response_load1.name_sji
index_el_sji_load1 = iris_response_load1.index_el_sji
geom_area_load1 = iris_response_load1.geom_area
elements_load1 = iris_response_load1.elements
comment_load1 = iris_response_load1.comment
version_load1 = iris_response_load1.version
date_load1 = iris_response_load1.date

# Version 2
# Tests for get_iris_response function
# Version 1
sav_file_path1 = os.path.join(rootdir, "idl_iris_get_response_20130903_new_version001.sav")
test_iris_response1 = scipy.io.readsav(sav_file_path1, python_dict=True, verbose=False)
iris_response_load1 = test_iris_response1["iris_response"][0]

lamb_load1 = iris_response_load1.lambda_vars
area_sg_load1 = iris_response_load1.area_sg
name_sg_load1 = iris_response_load1.name_sg
index_el_sg_load1 = iris_response_load1.index_el_sg
area_sji_load1 = iris_response_load1.area_sji
name_sji_load1 = iris_response_load1.name_sji
index_el_sji_load1 = iris_response_load1.index_el_sji
geom_area_load1 = iris_response_load1.geom_area
elements_load1 = iris_response_load1.elements
comment_load1 = iris_response_load1.comment
version_load1 = iris_response_load1.version
date_load1 = iris_response_load1.date

sav_file_path2 = os.path.join(rootdir, "idl_iris_get_response_20130903_new_version002.sav")
test_iris_response2 = scipy.io.readsav(sav_file_path2, python_dict=True, verbose=False)
iris_response_load2 = test_iris_response2["iris_response"][0]

lamb_load2 = iris_response_load2.lambda_vars
area_sg_load2 = iris_response_load2.area_sg
name_sg_load2 = iris_response_load2.name_sg
index_el_sg_load2 = iris_response_load2.index_el_sg
area_sji_load2 = iris_response_load2.area_sji
name_sji_load2 = iris_response_load2.name_sji
index_el_sji_load2 = iris_response_load2.index_el_sji
elements_load2 = iris_response_load2.elements
comment_load2 = iris_response_load2.comment
version_load2 = iris_response_load2.version
date_load2 = iris_response_load2.date

# Version 3
sav_file_path3 = os.path.join(rootdir, "idl_iris_get_response_20130903_new_version003.sav")

test_iris_response3 = scipy.io.readsav(sav_file_path3, python_dict=True, verbose=False)
iris_response_load3 = test_iris_response3["iris_response"][0]

date_obs_load3 = iris_response_load3.date_obs
lamb_load3 = iris_response_load3.lambda_vars
area_sg_load3 = iris_response_load3.area_sg
name_sg_load3 = iris_response_load3.name_sg
dn2phot_sg_load3 = iris_response_load3.dn2phot_sg
area_sji_load3 = iris_response_load3.area_sji
name_sji_load3 = iris_response_load3.name_sji
dn2phot_sji_load3 = iris_response_load3.dn2phot_sji
comment_load3 = iris_response_load3.comment
version_load3 = iris_response_load3.version
version_date_load3 = iris_response_load3.version_date

# Version 4
sav_file_path4 = os.path.join(rootdir, "idl_iris_get_response_20130903_new_version004.sav")
test_iris_response4 = scipy.io.readsav(sav_file_path4, python_dict=True, verbose=False)
iris_response_load4 = test_iris_response4["iris_response"][0]

date_obs_load4 = iris_response_load4.date_obs
lamb_load4 = iris_response_load4.lambda_vars
area_sg_load4 = iris_response_load4.area_sg
name_sg_load4 = iris_response_load4.name_sg
dn2phot_sg_load4 = iris_response_load4.dn2phot_sg
area_sji_load4 = iris_response_load4.area_sji
name_sji_load4 = iris_response_load4.name_sji
dn2phot_sji_load4 = iris_response_load4.dn2phot_sji
comment_load4 = iris_response_load4.comment
version_load4 = iris_response_load4.version
version_date_load4 = iris_response_load4.version_date


@pytest.fixture
@pytest.mark.remote_data
def iris_response1():
    # For testing of version 1
    return utils.get_iris_response(time_obs=parse_time("2013-09-03", format="utime"), response_version=1)


@pytest.mark.remote_data
def test_get_iris_response_version1(iris_response1):
    np_test.assert_almost_equal(iris_response1["AREA_SG"].value, area_sg_load1, decimal=6)
    np_test.assert_almost_equal(iris_response1["AREA_SJI"].value, area_sji_load1, decimal=6)


@pytest.fixture
@pytest.mark.remote_data
def iris_response2():
    # For testing of version 2
    return utils.get_iris_response(time_obs=parse_time("2013-09-03", format="utime"), response_version=2)


@pytest.mark.remote_data
def test_get_iris_response_version2(iris_response2):
    np_test.assert_almost_equal(iris_response2["AREA_SG"].value, area_sg_load2, decimal=6)
    np_test.assert_almost_equal(iris_response2["AREA_SJI"].value, area_sji_load2, decimal=6)


@pytest.fixture
@pytest.mark.remote_data
def iris_response3():
    # For testing of version 3
    return utils.get_iris_response(time_obs=parse_time("2013-09-03", format="utime"), response_version=3)


@pytest.mark.remote_data
def test_get_iris_response_version3(iris_response3):
    np_test.assert_almost_equal(iris_response3["AREA_SG"].value, area_sg_load3, decimal=6)
    np_test.assert_almost_equal(iris_response3["AREA_SJI"].value, area_sji_load3, decimal=6)


@pytest.fixture
@pytest.mark.remote_data
def iris_response4():
    # For testing of version 4
    return utils.get_iris_response(time_obs=parse_time("2013-09-03", format="utime"), response_version=4)


@pytest.mark.remote_data
def test_get_iris_response_version4(iris_response4):
    np_test.assert_almost_equal(iris_response4["AREA_SG"].value, area_sg_load4, decimal=3)
    np_test.assert_almost_equal(iris_response4["AREA_SJI"].value, area_sji_load4, decimal=3)


@pytest.mark.parametrize(
    "input_arrays, expected_array",
    [
        ([time_obs.value, time_cal_coeffs0, cal_coeffs0], iris_fit_expected0),
        ([time_obs.value, time_cal_coeffs0, cal_coeffs1], iris_fit_expected1),
        ([time_obs.value, time_cal_coeffs0, cal_coeffs2], iris_fit_expected2),
        ([time_obs.value, time_cal_coeffs1, cal_coeffs3], iris_fit_expected3),
        ([time_obs.value, time_cal_coeffs1, cal_coeffs4], iris_fit_expected4),
        ([time_obs.value, time_cal_coeffs1, cal_coeffs5], iris_fit_expected5),
    ],
)
def test_fit_iris_xput(input_arrays, expected_array):
    np_test.assert_almost_equal(
        utils.fit_iris_xput(input_arrays[0], input_arrays[1], input_arrays[2]),
        expected_array,
        decimal=6,
    )
