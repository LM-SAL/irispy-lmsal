"""
==================
Mg II Dopplergrams
==================

In this example we are going to produce a Dopplergram for the Mg II k line from a IRIS 400-step raster.
The Dopplergram is obtained by subtracting the intensities at symmetrical velocity shifts from the line core
(e.g. ±50 km/s).
For this kind of analysis we need a consistent wavelength calibration for each step of the raster.

This very large dense raster took more than three hours to complete the 400 scans (30 s exposures), which
means that the spacecraft's orbital velocity changes during the observations.
This means that any precise wavelength calibration will need to correct for those shifts.
"""
# sphinx_gallery_thumbnail_number = 5

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SpectralCoord
from astropy.io import fits
from scipy.constants import c
from scipy.interpolate import interp1d

from irispy.io import read_spectrograph_lvl2
from irispy.utils import image_clipping
from irispy.utils.utils import _download_data

###############################################################################
# We start with getting the data. This is done by downloading the data from the IRIS archive.
#
# In this case, we will use requests as to keep this example self contained
# but using your browser will also work.
#
# Using the url: http://www.lmsal.com/solarsoft/irisa/data/level2_compressed/2014/07/08/
# we are after `iris_l2_20140708_114109_3824262996_raster.tar.gz` which is ~730 MB in size.


urls = [
    "http://www.lmsal.com/solarsoft/irisa/data/level2_compressed/2014/07/08/20140708_114109_3824262996/iris_l2_20140708_114109_3824262996_raster.tar.gz"
]
_download_data(urls)
raster_filename = "iris_l2_20140708_114109_3824262996_raster_t000_r00000.fits"

###############################################################################
# We will open file using `irispy`.
# Note that when ``memmap=True``, the data values are read from the FITS file
# directly without the scaling to Float32, the data values are no longer in DN,
# but in scaled integer units that start at −2$^{16}$/2.

raster = read_spectrograph_lvl2(raster_filename, memmap=True, uncertainty=False)
print(raster)

###############################################################################
# We see that Mg II k 2796 is the last key.
# We will plot the spatially averaged spectrum

mg_ii = raster["Mg II k 2796"][0]
(mg_wave,) = mg_ii.axis_world_coords("wl")

plt.plot(mg_wave.to("nm"), mg_ii.data.mean((0, 1)))

plt.show()

###############################################################################
# To better understand the orbital velocity problem let us look at how the
# line intensity varies for a strong Mn I line at around 280.2 nm, in
# between the Mg II k and h lines. For this dataset, the line core of this
# line falls around index 350. To plot a spectroheliogram in the correct
# orientation we will transpose the data

lower_corner = [SpectralCoord(280.2, unit=u.nm), None]
upper_corner = [SpectralCoord(280.2, unit=u.nm), None]
mg_crop = mg_ii.crop(lower_corner, upper_corner)

plt.figure(figsize=(6, 10))
vmin, vmax = image_clipping(mg_ii.data[..., 350])
ax = mg_crop.plot(vmin=vmin, vmax=vmax + 1000)

plt.show()

###############################################################################
# You can see a regular bright-dark pattern along the x axis, an
# indication that its intensities are not taken at the same position in
# the line because of wavelength shifts. The shifts are caused by the
# orbital velocity changes, and we can find these in the auxiliary
# metadata which are to be found in the extension past the "last" window
# in the data set

aux = fits.getdata(raster_filename, 9)
aux_hd = fits.getheader(raster_filename, 9)
v_obs = aux[:, aux_hd["OBS_VRIX"]]
v_obs /= 1000.0  # convert to km/s

plt.plot(v_obs)
plt.ylabel("Orbital velocity (km/s)")
plt.xlabel("Scan number")

plt.show()

###############################################################################
# To look at intensities at any given scan we only need to subtract this
# velocity shift from the wavelength scale, but to look at the whole image
# at a given wavelength we must interpolate the original data to take this
# shift into account. Here is a way to do it (note that array dimensions
# apply to this specific set only!)

c_kms = c / 1000.0
wave_shift = -v_obs * mg_wave[350] / (c_kms)
# Linear interpolation in wavelength, for each scan
for i in range(mg_ii.data.shape[0]):
    tmp = interp1d(mg_wave - wave_shift[i], mg_ii.data[:, i, :], bounds_error=False)
    mg_ii.data[:, i, :] = tmp(mg_wave)

###############################################################################
# Now we can plot the shifted data to see that the large scale shifts
# have disappeared

plt.figure(figsize=(6, 10))
vmin, vmax = image_clipping(mg_ii.data[..., 350])
plt.imshow(
    mg_ii.data[..., 350].T,
    origin="lower",
    aspect=0.5,
    vmin=vmin,
    vmax=vmax,
)

plt.show()

###############################################################################
# Some residual shift remains, but we will not correct for it here. A more
# elaborate correction can be obtained by the IDL routine
# ``iris_prep_wavecorr_l2``, but this has not yet been ported to Python
# see the `IDL version of this
# tutorial <http://iris.lmsal.com/itn26/tutorials.html#mg-ii-dopplergrams>`__
# for more details.
#
# We can use the calibrated data for example to calculate Dopplergrams. A
# Dopplergram is here defined as the difference between the intensities at
# two wavelength positions at the same (and opposite) distance from the
# line core. For example, at +/- 50 km/s from the Mg II k3 core. To do
# this, let us first calculate a velocity scale for the k line and find
# the indices of the -50 and +50 km/s velocity positions (here using the
# convention of negative velocities for up flows):

mg_k_centre = 279.6351 * u.nm
pos = 50  # in km/s around line centre
velocity = (mg_wave - mg_k_centre) * c_kms / mg_k_centre
index_p = np.argmin(np.abs(velocity - pos))
index_m = np.argmin(np.abs(velocity + pos))
doppl = mg_ii.data[..., index_m] - mg_ii.data[..., index_p]

###############################################################################
# And now we can plot this as before (intensity units are again arbitrary
# because of the unscaled DNs):

fig = plt.figure(figsize=(6, 10))
vmin, vmax = image_clipping(doppl)
plt.imshow(
    doppl.T,
    cmap="gist_gray",
    origin="lower",
    aspect=0.5,
    vmin=vmin,
    vmax=vmax,
)

plt.show()
