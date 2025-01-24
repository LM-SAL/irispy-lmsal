"""
================
Spectral fitting
================

In this example, we are going to fit spectral lines from IRIS data with Gaussians.
We use this to work out the moments.
"""

import matplotlib.pyplot as plt
import numpy as np
import pooch

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.modeling import models as m
from astropy.modeling.fitting import TRFLSQFitter, parallel_fit_dask
from astropy.visualization import time_support

from sunpy.coordinates.frames import Helioprojective

from irispy.io import read_files

time_support()

###############################################################################
# We will start by getting some data from the IRIS archive.
#
# In this case, we will use ``pooch`` so to keep this example self-contained
# but using your browser will also work.
#
# Using the url: http://www.lmsal.com/solarsoft/irisa/data/level2_compressed/2013/09/02/20130902_163935_4000255147/
# we are after the raster sequence (~900 MB).

raster_filename = pooch.retrieve(
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

raster = read_files(raster_filename, memmap=False)

###############################################################################
# We will just focus on the Si IV 1403 line which we can select using a key.
# Then we will just plot a spectral line selected at random in space.

# There is only one complete scan, so we index that away.
si_iv_1403 = raster["Si IV 1403"][0]

lower_corner = [None, SkyCoord(-338 * u.arcsec, 275 * u.arcsec, frame=Helioprojective)]
upper_corner = [None, SkyCoord(-338 * u.arcsec, 275 * u.arcsec, frame=Helioprojective)]

si_iv_crop = si_iv_1403.crop(lower_corner, upper_corner)

fig = plt.figure()
ax = fig.add_subplot(111, projection=si_iv_crop.wcs)
si_iv_crop.plot(axes=ax)

###############################################################################
# We will want to make two rebinned cubes from the full raster,
# one summed along the wavelength dimension and one of the spectra averaged over all spatial pixels.

wl_sum = si_iv_1403.rebin((1, 1, si_iv_1403.data.shape[-1]), operation=np.sum)[0]
print(wl_sum)

spatial_mean = si_iv_1403.rebin((*si_iv_1403.data.shape[:-1], 1))[0, 0, :]
print(spatial_mean)

# ###############################################################################
# So this is what the spatial mean looks like.

fig = plt.figure()
ax = fig.add_subplot(111, projection=spatial_mean.wcs)
ax = spatial_mean.plot(axes=ax)

# ###############################################################################
# Now we can create a model for this spectra.

si_iv_core = 140.28 * u.nm

initial_model = m.Const1D(amplitude=2 * si_iv_1403.unit) + m.Gaussian1D(
    amplitude=8 * si_iv_1403.unit, mean=si_iv_core, stddev=0.005 * u.nm
)
print(initial_model)

# ###############################################################################
# To improve our initial conditions we now fit the initial model to the spatially averaged spectra.
# To do this we use the `.axis_world_coords` method of NDCube which returns all,
# or a subset of the world coordinates along however many array axes they are
# correlated with. So in this case we get the wavelength dimension which only
# returns a single `SpectralCoord` object corresponding to the first array dimension of the cube.

fitter = TRFLSQFitter(calc_uncertainties=True)
average_fit = fitter(
    initial_model,
    spatial_mean.axis_world_coords("em.wl")[0].to(u.nm),
    spatial_mean.data * spatial_mean.unit,
)
print(average_fit)

# ###############################################################################
# Now we can add to our previous plot, the initial model and the model fit to the average spectra.

fig = plt.figure()
ax = spatial_mean.plot(label="Spatial average")
ax.plot(initial_model(spatial_mean.axis_world_coords("em.wl")[0].to(u.nm)), label="Initial model")
ax.plot(average_fit(spatial_mean.axis_world_coords("em.wl")[0].to(u.nm)), linestyle="--", label="Spatial average fit")
plt.legend()

# ###############################################################################
# Now we have our model to fit to all of the spectra, for all slit steps.
#
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

# We can therefore fit the cube as follows:
spice_model_fit = parallel_fit_dask(
    data=si_iv_1403[0],
    data_unit=si_iv_1403.unit,
    model=average_fit,
    fitter=TRFLSQFitter(),
    fitting_axes=1,
    # Filter out non-finite values otherwise the fitting will fail
    fitter_kwargs={"filter_non_finite": True},
)

# # ###############################################################################
# # Given that we are going to want to visualise the output of a few fits, I am going to define a plotting function which will display the shift in the peak locations of the two Gaussians. We shall talk more about this later.


def plot_spice_fit(spice_model_fit):
    g1_peak_shift = spice_model_fit.mean_1.quantity.to(u.km / u.s, equivalencies=u.doppler_optical(NIV_wave))
    g2_peak_shift = spice_model_fit.mean_2.quantity.to(u.km / u.s, equivalencies=u.doppler_optical(NeVIII_wave))

    fig, axs = plt.subplots(ncols=3, subplot_kw={"projection": wl_sum}, figsize=(11, 4))
    fig.suptitle(f"SPICE - {hdu.header['EXTNAME']} - {hdu.header['DATE-AVG']}")

    wl_sum.plot(axes=axs[0])
    fig.colorbar(axs[0].get_images()[0], ax=axs[0], extend="both", label=f"{wl_sum.unit:latex}", shrink=0.8)
    axs[0].set_title("Data (summed over wavelength)", pad=40)

    g1_max = np.percentile(np.abs(g1_peak_shift.value), 99)
    mean_1 = axs[1].imshow(g1_peak_shift.value, cmap="coolwarm", vmin=-g1_max, vmax=g1_max)
    fig.colorbar(
        mean_1, ax=axs[1], extend="both", label=f"Velocity from Doppler shift [{g1_peak_shift.unit:latex}]", shrink=0.8
    )
    axs[1].set_title(f"N IV ({NIV_wave:latex})", pad=40)

    g2_max = np.percentile(np.abs(g2_peak_shift.value), 98)
    mean_2 = axs[2].imshow(g2_peak_shift.value, cmap="coolwarm", vmin=-g2_max, vmax=g2_max)
    fig.colorbar(
        mean_2, ax=axs[2], extend="both", label=f"Velocity from Doppler shift [{g2_peak_shift.unit:latex}]", shrink=0.8
    )
    axs[2].set_title(f"Ne VIII ({NeVIII_wave:latex})", pad=40)

    for ax in axs:
        ax.coords[0].set_ticklabel(exclude_overlapping=True)
        ax.coords[0].set_axislabel("Helioprojective Longitude")
        ax.coords[1].set_axislabel("Helioprojective Latitude")
        ax.coords[2].set_ticklabel(exclude_overlapping=True)

    fig.tight_layout()


plot_spice_fit(spice_model_fit)


plt.show()
