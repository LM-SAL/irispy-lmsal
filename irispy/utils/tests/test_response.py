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
time_obs = [parse_time("2013-09-03")]
time_cal_coeffs0 = parse_time(iris_response.get("C_F_TIME"), format="utime")
time_cal_coeffs1 = parse_time(iris_response.get("C_N_TIME"), format="utime")
cal_coeffs0 = iris_response.get("COEFFS_FUV")[0, :, :]
cal_coeffs1 = iris_response.get("COEFFS_FUV")[1, :, :]
cal_coeffs2 = iris_response.get("COEFFS_FUV")[2, :, :]
cal_coeffs3 = iris_response.get("COEFFS_NUV")[0, :, :]
cal_coeffs4 = iris_response.get("COEFFS_NUV")[1, :, :]
cal_coeffs5 = iris_response.get("COEFFS_NUV")[2, :, :]
cal_coeffs6 = iris_response.get("COEFFS_NUV")[3, :, :]
iris_fit_expected0 = np.array([0.79865301])
iris_fit_expected1 = np.array([2.2495413])
iris_fit_expected2 = np.array([2.2495413])
iris_fit_expected3 = np.array([0.23529011])
iris_fit_expected4 = np.array([0.25203046])
iris_fit_expected5 = np.array([0.25265095])
iris_fit_expected6 = np.array([0.253833])


def test_get_iris_response_greater_than_6():
    with pytest.raises(KeyError):
        utils.get_iris_response(time_obs, response_version=13)


# Tests for get_iris_response function
# Version 1
sav_file_path1 = os.path.join(rootdir, "idl_iris_get_response_20130903_new_version001.sav")
test_iris_response1 = scipy.io.readsav(sav_file_path1, python_dict=True, verbose=False)
iris_response_load1 = test_iris_response1["iris_response"][0]
area_sg_load1 = iris_response_load1.area_sg
area_sji_load1 = iris_response_load1.area_sji

# Version 2
sav_file_path2 = os.path.join(rootdir, "idl_iris_get_response_20130903_new_version002.sav")
test_iris_response2 = scipy.io.readsav(sav_file_path2, python_dict=True, verbose=False)
iris_response_load2 = test_iris_response2["iris_response"][0]
area_sg_load2 = iris_response_load2.area_sg
area_sji_load2 = iris_response_load2.area_sji

# Version 3
sav_file_path3 = os.path.join(rootdir, "idl_iris_get_response_20130903_new_version003.sav")

test_iris_response3 = scipy.io.readsav(sav_file_path3, python_dict=True, verbose=False)
iris_response_load3 = test_iris_response3["iris_response"][0]
area_sg_load3 = iris_response_load3.area_sg
area_sji_load3 = iris_response_load3.area_sji

# Version 4
sav_file_path4 = os.path.join(rootdir, "idl_iris_get_response_20130903_new_version004.sav")
test_iris_response4 = scipy.io.readsav(sav_file_path4, python_dict=True, verbose=False)
iris_response_load4 = test_iris_response4["iris_response"][0]
area_sg_load4 = iris_response_load4.area_sg
area_sji_load4 = iris_response_load4.area_sji

# Version 5
sav_file_path5 = os.path.join(rootdir, "idl_iris_get_response_20130903_new_version005.sav")
test_iris_response5 = scipy.io.readsav(sav_file_path5, python_dict=True, verbose=False)
iris_response_load5 = test_iris_response5["iris_response"][0]
area_sg_load5 = iris_response_load5.area_sg
area_sji_load5 = iris_response_load5.area_sji

# Version 6
sav_file_path6 = os.path.join(rootdir, "idl_iris_get_response_20130903_new_version006.sav")
test_iris_response6 = scipy.io.readsav(sav_file_path6, python_dict=True, verbose=False)
iris_response_load6 = test_iris_response6["iris_response"][0]
area_sg_load6 = iris_response_load6.area_sg
area_sji_load6 = iris_response_load6.area_sji


@pytest.fixture
def iris_responsev1():
    return utils.get_iris_response(time_obs=parse_time("2013-09-03"), response_version=1)


def test_get_iris_response_version1(iris_responsev1):
    np_test.assert_almost_equal(iris_responsev1["AREA_SG"].value, area_sg_load1, decimal=6)
    np_test.assert_almost_equal(iris_responsev1["AREA_SJI"].value, area_sji_load1, decimal=6)


@pytest.fixture
def iris_responsev2():
    return utils.get_iris_response(time_obs=parse_time("2013-09-03"), response_version=2)


def test_get_iris_response_version2(iris_responsev2):
    np_test.assert_almost_equal(iris_responsev2["AREA_SG"].value, area_sg_load2, decimal=6)
    np_test.assert_almost_equal(iris_responsev2["AREA_SJI"].value, area_sji_load2, decimal=6)


@pytest.fixture
def iris_responsev3():
    return utils.get_iris_response(time_obs=parse_time("2013-09-03"), response_version=3)


def test_get_iris_response_version3(iris_responsev3):
    np_test.assert_almost_equal(iris_responsev3["AREA_SG"].value, area_sg_load3, decimal=6)
    np_test.assert_almost_equal(iris_responsev3["AREA_SJI"].value, area_sji_load3, decimal=6)


@pytest.fixture
def iris_responsev4():
    return utils.get_iris_response(time_obs=parse_time("2013-09-03"), response_version=4)


@pytest.fixture
def iris_responsev5():
    return utils.get_iris_response(time_obs=parse_time("2013-09-03"), response_version=5)


@pytest.fixture
def iris_responsev6():
    return utils.get_iris_response(time_obs=parse_time("2013-09-03"), response_version=6)


def test_get_iris_response_version4(iris_responsev4):
    np_test.assert_almost_equal(iris_responsev4["AREA_SG"].value, area_sg_load4, decimal=3)
    np_test.assert_almost_equal(iris_responsev4["AREA_SJI"].value, area_sji_load4, decimal=3)


def test_get_iris_response_version5(iris_responsev5):
    np_test.assert_almost_equal(iris_responsev5["AREA_SG"].value, area_sg_load5, decimal=3)
    np_test.assert_almost_equal(iris_responsev5["AREA_SJI"].value, area_sji_load5, decimal=3)


def test_get_iris_response_version6(iris_responsev6):
    np_test.assert_almost_equal(iris_responsev6["AREA_SG"].value, area_sg_load6, decimal=3)
    np_test.assert_almost_equal(iris_responsev6["AREA_SJI"].value, area_sji_load6, decimal=3)


@pytest.mark.parametrize(
    "input_arrays, expected_array",
    [
        ([time_obs, time_cal_coeffs0, cal_coeffs0], iris_fit_expected0),
        ([time_obs, time_cal_coeffs0, cal_coeffs1], iris_fit_expected1),
        ([time_obs, time_cal_coeffs0, cal_coeffs2], iris_fit_expected2),
        ([time_obs, time_cal_coeffs1, cal_coeffs3], iris_fit_expected3),
        ([time_obs, time_cal_coeffs1, cal_coeffs4], iris_fit_expected4),
        ([time_obs, time_cal_coeffs1, cal_coeffs5], iris_fit_expected5),
        ([time_obs, time_cal_coeffs1, cal_coeffs6], iris_fit_expected6),
    ],
)
def test_fit_iris_xput(input_arrays, expected_array):
    np_test.assert_almost_equal(
        utils.fit_iris_xput(input_arrays[0], input_arrays[1], input_arrays[2]),
        expected_array,
        decimal=6,
    )


def test_fixed_time():
    sav_file_path6 = os.path.join(rootdir, "test_response.sav")
    response = scipy.io.readsav(sav_file_path6, python_dict=True, verbose=False)
    response = response["test"][0]
    area_sg = response.area_sg
    area_sji = response.area_sji
    python_response = utils.get_iris_response(parse_time("2013-08-31"), response_version=6)
    np_test.assert_almost_equal(python_response["AREA_SG"].value, area_sg, decimal=3)
    np_test.assert_almost_equal(python_response["AREA_SJI"].value, area_sji, decimal=3)
