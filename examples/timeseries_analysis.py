"""
====================
Time series analysis
====================

In this example, we are going to work with spectra and slit-jaw images to study dynamical phenomena.
The subject of this example is umbral flashes.
"""
# sphinx_gallery_thumbnail_number = 4


import astropy.units as u
import matplotlib.pyplot as plt
from astropy.coordinates import SpectralCoord
from astropy.visualization import time_support

from irispy.io import read_sji_lvl2, read_spectrograph_lvl2
from irispy.utils import image_clipping
from irispy.utils.utils import _download_data

time_support()

###############################################################################
# We start with getting the data. This is done by downloading the data from the IRIS archive.
#
# In this case, we will use requests as to keep this example self contained
# but using your browser will also work.
#
# Using the url: http://www.lmsal.com/solarsoft/irisa/data/level2_compressed/2013/09/02/20130902_163935_4000255147/
# we are after the 1400 slit-jaw file and the raster file (~900 MB).

urls = [
    "http://www.lmsal.com/solarsoft/irisa/data/level2_compressed/2013/09/02/20130902_163935_4000255147/iris_l2_20130902_163935_4000255147_SJI_1400_t000.fits.gz",
    "http://www.lmsal.com/solarsoft/irisa/data/level2_compressed/2013/09/02/20130902_163935_4000255147/iris_l2_20130902_163935_4000255147_raster.tar.gz",
]
_download_data(urls)
raster_filename = "iris_l2_20130902_163935_4000255147_raster_t000_r00000.fits"
sji_filename = "iris_l2_20130902_163935_4000255147_SJI_1400_t000.fits.gz"

###############################################################################
# We will open file using `irispy`.
# Note that when ``memmap=True``, the data values are read from the FITS file
# directly without the scaling to Float32, the data values are no longer in DN,
# but in scaled integer units that start at −2$^{16}$/2.

raster = read_spectrograph_lvl2(raster_filename, memmap=True, uncertainty=False)
sji_1400 = read_sji_lvl2(sji_filename, memmap=True, uncertainty=False)
print(raster)

###############################################################################
# We see that Mg II k 2796 is the last key.
# For this dataset the spectral cadence is about 3 seconds. The Mg II k3
# core is located around wavelength pixel 103. We can use this information
# to make a space-time image of the Mg II k3 wavelength.

mg_ii = raster["Mg II k 2796"][0]
c_ii = raster["C II 1336"][0]

lower_corner = [SpectralCoord(279.63, unit=u.nm), None]
upper_corner = [SpectralCoord(279.63, unit=u.nm), None]
mg_crop = mg_ii.crop(lower_corner, upper_corner)

vmin, vmax = image_clipping(mg_ii.data[..., 103])
mg_crop.plot(vmin=vmin, vmax=vmax)
plt.xlabel("Solar Y")

plt.show()

###############################################################################
# The middle section between 60"-75" is on the umbra of a sunspot, even though
# it is not obvious from this image. One can see very clearly the umbral oscillations,
# with a clear regular pattern of dark/bright streaks.
#
# Let us now load the 1400 slit-jaw and plot it for context

vmin, vmax = image_clipping(sji_1400.data[0])
sji_1400.plot(vmin=vmin, vmax=vmax)

plt.show()

###############################################################################
# The slit pixel 220 is a location on the sunspot's umbra. We will use it
# to get some plots. For example, let's plot the k3 intensity (spectral
# pixel 103 of ``mg_ii``) and the core of the brightest C II line
# (spectral pixel 90 of ``c_ii``) vs. time in minutes (showing first ~10
# minutes only)

mg_ii_times = mg_ii.axis_world_coords("time", wcs=mg_ii.extra_coords)[0][:200]
c_ii_times = c_ii.axis_world_coords("time", wcs=c_ii.extra_coords)[0][:200]
plt.plot(mg_ii_times, mg_ii.data[:200, 220, 103])
plt.plot(c_ii_times, c_ii.data[:200, 220, 90])

plt.show()

###############################################################################
# Imagine now you wanted to compare these oscillations with
# the intensity on the slit-jaw. How to do it? The slit-jaws are typically
# taken at a different cadence, so you will need to load the corresponding
# time array for the 1400 slit-jaw

times_sji = sji_1400.axis_world_coords("time", wcs=sji_1400.extra_coords)[0][:50]

###############################################################################
# Now we can plot both spectral lines and slit-jaw for a pixel close
# to the slit at the same y position (index 220)

plt.plot(mg_ii_times, mg_ii.data[:200, 220, 103])
plt.plot(c_ii_times, c_ii.data[:200, 220, 90])
plt.plot(times_sji, sji_1400.data[:50, 190, 220], "y-")

plt.show()

###############################################################################
# You are now ready to explore all the correlations, anti-correlations,
# and phase differences.
