"""
================
Spectral fitting
================

In this example, we are going to fit spectral lines from IRIS, using the raster data with a single Gaussian.
Then use the fitted values to calculate the Gaussian moments.

This example is laid out like this so it can be run as a script without erroring.
If you want to run it in a notebook, you can just copy the code below
and paste it into a cell, you do not need the functions or the ``if __name__ == "__main__":`` part.
"""

import shutil
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pooch
from dask.distributed import Client

import astropy.units as u
from astropy import constants
from astropy.coordinates import SkyCoord, SpectralCoord
from astropy.modeling import models as m
from astropy.modeling.fitting import LevMarLSQFitter, parallel_fit_dask
from astropy.visualization import time_support

from sunpy.coordinates import frames

from irispy.io import read_files

###############################################################################
# We will start by getting some data from the IRIS archive.
#
# In this case, we will use ``pooch`` so to keep this example self-contained
# but using your browser will also work.
#
# Using the url: http://www.lmsal.com/solarsoft/irisa/data/level2_compressed/2018/01/02/20180102_153155_3610108077/iris_l2_20180102_153155_3610108077_raster.tar.gz
# we are after the raster sequence (~300 MB).


def get_data():
    return pooch.retrieve(
        "http://www.lmsal.com/solarsoft/irisa/data/level2_compressed/2018/01/02/20180102_153155_3610108077/iris_l2_20180102_153155_3610108077_raster.tar.gz",
        known_hash="0ec2b7b20757c52b02e0d92c27a5852b6e28759512c3d455f8b6505d4e1f5cd6",
    )


###############################################################################
# Now to open the files using ``irispy-lmsal``.

# Note that when ``memmap=True``, the data values are read from the FITS file
# directly without the scaling to Float32 (via "b_zero" and "b_scale"),
# the data values are no longer in DN, but in scaled integer units that start at -2$^{16}$/2.
#
# We will use ``memmap=False`` because we want to fit the actual the data values.


def read_file(raster_filename):
    return read_files(raster_filename, memmap=False)


###############################################################################
# We will just focus on the Si IV 1403 line which we can select using a key.
# Then we will just plot a spectral line selected at random in space.


def get_si_iv_1403(raster):
    # There is only one complete scan, so we index that away.
    si_iv_1403 = raster["Si IV 1403"][0]

    # However, before we get to that, we will shrink the data cube to make it easier to work with.
    top_left = [None, SkyCoord(-290 * u.arcsec, 260 * u.arcsec, frame=frames.Helioprojective)]
    bottom_right = [None, SkyCoord(-360 * u.arcsec, 310 * u.arcsec, frame=frames.Helioprojective)]
    return si_iv_1403.crop(top_left, bottom_right)


###############################################################################
# Let us just check the full field of view at the line core.


def quicklook_si_iv_1403(si_iv_1403, si_iv_core):
    lower_corner = [SpectralCoord(si_iv_core), None]
    upper_corner = [SpectralCoord(si_iv_core), None]
    si_iv_spec_crop = si_iv_1403.crop(lower_corner, upper_corner)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection=si_iv_spec_crop.wcs)
    ax = si_iv_spec_crop.plot(axes=ax, vmin=0, vmax=100)

    return si_iv_spec_crop, fig


###############################################################################
# We will want to make two rebinned cubes from the full raster,
# one summed along the wavelength dimension and one of the spectra averaged over all spatial pixels.


def get_averages(si_iv_1403):
    wl_sum = si_iv_1403.rebin((1, 1, si_iv_1403.data.shape[-1]), operation=np.sum)[0]
    spatial_mean = si_iv_1403.rebin((*si_iv_1403.data.shape[:-1], 1))[0, 0, :]
    return wl_sum, spatial_mean


################################################################################
# Now we can create a model for this spectra.


def setup_initial_model(si_iv_1403, si_iv_core):
    return m.Const1D(amplitude=2 * si_iv_1403.unit) + m.Gaussian1D(
        amplitude=8 * si_iv_1403.unit, mean=si_iv_core, stddev=0.005 * u.nm
    )


################################################################################
# To improve our initial conditions we now fit the initial model to the spatially averaged spectra.
# To do this we use the `ndcube.NDCube.axis_world_coords` method of NDCube which returns all,
# or a subset of the world coordinates along however many array axes they are
# correlated with. So in this case we get the wavelength dimension which only
# returns a single `SpectralCoord` object corresponding to the first array dimension of the cube.


def get_average_fit(initial_model, spatial_mean):
    fitter = LevMarLSQFitter(calc_uncertainties=True)
    return fitter(
        initial_model,
        spatial_mean.axis_world_coords("em.wl")[0].to(u.nm),
        spatial_mean.data * spatial_mean.unit,
    )


################################################################################
# Now we check, the initial model and the model fitted to the average spectra.


def plot_initial_model(initial_model, spatial_mean, average_fit):
    fig = plt.figure()
    ax = spatial_mean.plot(label="Spatial average")
    ax.plot(initial_model(spatial_mean.axis_world_coords("em.wl")[0].to(u.nm)), label="Initial model")
    ax.plot(
        average_fit(spatial_mean.axis_world_coords("em.wl")[0].to(u.nm)), linestyle="--", label="Spatial average fit"
    )
    plt.legend()
    return fig


################################################################################
# The function `.parallel_fit_dask` will map a model to each element of a cube along
# one (or more) "fitting axes", in this case our fitting axis is our wavelength
# axis (array axis -1). So we want to fit each slice of the data array along the 3rd axis.
#
# The key arguments to the parallel_fit_dask function are:
#
# * A data array: This can be a numpy array or a dask array, or a NDData (or subclass like NDCube)
#                 object. If it's a NDData object then the data, wcs, mask, data_unit and uncertainty
#                 are all extracted from the NDData object and used in place of their respective keyword
#                 arguments.
# * A model to fit
# * A fitter instance.
# * The fitting axis (or axes).
#
# What is returned from `.parallel_fit_dask` is a model with array parameters with
# the shape of the non-fitting axes of the data.


def do_fitting(si_iv_1403, average_fit):
    """
    This function will do the fitting.

    It is not needed if you are running this in a notebook.
    """
    # First we will start a local dask client to improve the fitting performance
    client = Client()

    # Next we will want to setup error reporting
    # Since it is using multiprocess/dask, logs are returned instead of errors
    diag_path = Path("./diag")
    shutil.rmtree(diag_path, ignore_errors=True)

    # LevMarLSQFitter does not like negative values, so we need to filter the data
    # You can also set a mask but for now, we will just handle the data.
    filtered_data = np.where(si_iv_1403.data < 0, 0, si_iv_1403.data)
    # We also need to set any NaN nor non-finite values to 0, otherwise the fitting will fail
    filtered_data = np.where(~np.isfinite(filtered_data), 0, filtered_data)

    # We can therefore fit the cube as follows
    with warnings.catch_warnings():
        # There are several WCS warnings we just want to ignore
        warnings.simplefilter("ignore")
        iris_model_fit = parallel_fit_dask(
            data=filtered_data,
            data_unit=si_iv_1403.unit,
            fitting_axes=2,
            world=si_iv_1403.wcs,
            model=average_fit,
            fitter=LevMarLSQFitter(),
            diagnostics="error",
            diagnostics_path=diag_path,
            scheduler=client,
        )
    diag_folders = list(diag_path.glob("*"))
    errors = []
    for diag in diag_folders:
        if (path := (diag / "error.log")).exists():
            content = open(path).read()
            errors.append(content)
    print(f"{len(errors)} errors occurred")
    if errors:
        print("First error is:")
        print(errors[0])

    client.close()
    return iris_model_fit


################################################################################
# Given that we are going to want to visualize the output, a simple plotting
# function below has been created for this purpose.


def plot_gaussian_fit(model_fit, si_iv_spec_crop, si_iv_core):
    fig, axs = plt.subplots(nrows=1, ncols=3, subplot_kw={"projection": si_iv_spec_crop}, figsize=(16, 6))

    net_flux = (
        np.sqrt(2 * np.pi)
        * (model_fit.amplitude_0 + model_fit.amplitude_1)
        * model_fit.stddev_1.quantity
        / np.mean(si_iv_1403.axis_world_coords("wl")[0][1:] - si_iv_1403.axis_world_coords("wl")[0][:-1])
    )
    amp_max = np.nanpercentile(np.abs(net_flux.value), 99)
    amp = axs[0].imshow(net_flux.value, vmin=0, vmax=amp_max, origin="lower")
    cbar = fig.colorbar(amp, ax=axs[0])
    cbar.set_label(label=f"Intensity [{net_flux.unit.to_string()}]", fontsize=8)
    cbar.ax.tick_params(labelsize=8)
    axs[0].set_title("Gaussian Net Flux")

    core_shift = ((model_fit.mean_1.quantity.to(u.nm)) - si_iv_core) / si_iv_core * (constants.c.to(u.km / u.s))
    shift_max = np.nanpercentile(np.abs(core_shift.value), 95)
    shift = axs[1].imshow(core_shift.value, cmap="coolwarm", vmin=-shift_max, vmax=shift_max)
    cbar = fig.colorbar(shift, ax=axs[1], extend="both")
    cbar.set_label(label=f"Doppler shift [{core_shift.unit.to_string()}]", fontsize=8)
    cbar.ax.tick_params(labelsize=8)
    axs[1].set_title("Velocity from Gaussian shift")

    sigma = (model_fit.stddev_1.quantity.to(u.nm)) / si_iv_core * (constants.c.to(u.km / u.s))
    line_max = np.nanpercentile(np.abs(sigma.value), 95)
    line = axs[2].imshow(sigma.value, vmax=line_max)
    cbar = fig.colorbar(line, ax=axs[2])
    cbar.set_label(label=f"Line Width [{sigma.unit.to_string()}]", fontsize=8)
    cbar.ax.tick_params(labelsize=8)
    axs[2].set_title("Gaussian Sigma")

    for ax in axs:
        ax.coords[0].set_ticklabel(exclude_overlapping=True, fontsize=8)
        ax.coords[0].set_axislabel("Helioprojective Longitude", fontsize=8)
        ax.coords[1].set_ticklabel(exclude_overlapping=True, fontsize=8)
        ax.coords[1].set_axislabel("Helioprojective Latitude", fontsize=8)
    fig.tight_layout()
    return fig


################################################################################
# Let us see the output!

if __name__ == "__main__":
    time_support()

    raster_filename = get_data()
    raster = read_file(raster_filename)
    si_iv_1403 = get_si_iv_1403(raster)
    si_iv_core = 140.277 * u.nm
    si_iv_spec_crop, quicklook_fig = quicklook_si_iv_1403(si_iv_1403, si_iv_core)

    wl_sum, spatial_mean = get_averages(si_iv_1403)
    initial_model = setup_initial_model(si_iv_1403, si_iv_core)
    average_fit = get_average_fit(initial_model, spatial_mean)
    initial_fig = plot_initial_model(initial_model, spatial_mean, average_fit)

    iris_model_fit = do_fitting(si_iv_1403, average_fit)
    guass_fig = plot_gaussian_fit(iris_model_fit, si_iv_spec_crop, si_iv_core)
    plt.show()
