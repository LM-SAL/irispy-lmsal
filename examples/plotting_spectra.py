"""
================
Plotting Spectra
================

In this example, we will show how to plot the spectra of a given observation.

You can get IRIS data with co-aligned SDO data (and more) from https://iris.lmsal.com/search/
"""
# sphinx_gallery_thumbnail_number = 5

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import pooch
from astropy.coordinates import SkyCoord, SpectralCoord
from sunpy.coordinates.frames import Helioprojective

from irispy.io import read_files

###############################################################################
# We start with getting the data.
# This is done by downloading the data from the IRIS archive.
#
# In this case, we will use ``pooch`` as to keep this example self contained
# but using your browser will also work.

raster_filename = pooch.retrieve(
    "http://www.lmsal.com/solarsoft/irisa/data/level2_compressed/2018/01/02/20180102_153155_3610108077/iris_l2_20180102_153155_3610108077_raster.tar.gz",
    known_hash=None,
)

###############################################################################
# We will now open the raster file we just downloaded.

# Note that when ``memmap=True``, the data values are read from the FITS file
# directly without the scaling to Float32, the data values are no longer in DN,
# but in scaled integer units that start at −2$^{16}$/2.
raster = read_files(raster_filename, memmap=True)

# Printing will give us an overview of the file.
print(raster)
# We can get every wavelength window.
print(raster.keys())

###############################################################################
# We will now access one wavelength window and plot it.
# You can access it like indexing a dictionary.

# The index is needed to escape the data container.
mg_ii = raster["Mg II k 2796"][0]
# Overview of this window.
print(mg_ii)

###############################################################################
# Now, we will plot the full raster.

mg_ii.plot()

plt.show()

###############################################################################
# We can then slice again to fetch a specific time during the observation.
# We can index to get the first index and plot it.

print(mg_ii[0])

mg_ii[0].plot()

plt.show()

###############################################################################
# Or we can get the spectrum at a specific pixel.

mg_ii[120, 200].plot()

plt.show()

###############################################################################
# The default plot uses the units and labels from the WCS,
# and in this case, wavelength is in metres.
# To change this, we will do the following:

# This retrieves the full wavelength range.
(mg_wave,) = mg_ii.axis_world_coords("wl")
# We can change the units to what we need.
mg_wave.to("nm")

fig, ax = plt.subplots()
# This is an example of indexing the data array directly.
plt.plot(mg_wave.to("nm"), mg_ii.data[120, 200])

plt.show()

###############################################################################
# When we use the underlying data directly, we lose all the metadata and WCS.
# Therefore, suggest not plotting like this.
# We will improve on the default plot by adjusting a few options.

mg_ii[0].plot(aspect="auto", cmap="irissjiNUV")
ax = plt.gca()
ax.coords[0].set_major_formatter("x.x")
ax.coords[0].set_format_unit(u.nm)
ax.set_xlabel("Wavelength (nm)")
ax.set_ylabel("Helioprojective latitude [Solar Y] (arcsec)")
# Remove longitude ticks
ax.coords[2].set_ticks([] * u.degree)

plt.show()

###############################################################################
# What is the wavelength position that corresponds to Mg II k core (279.63 nm)?

wcs_loc = mg_ii.wcs.world_to_pixel(
    SpectralCoord(279.63, unit=u.nm),
    SkyCoord(0 * u.arcsec, 0 * u.arcsec, frame=Helioprojective),
)
mg_index = int(np.round(wcs_loc[0]))

print(mg_index)

###############################################################################
# Now we will plot spectra for Mg II k core wavelength.
# This uses the ``crop`` method.

lower_corner = [SpectralCoord(280, unit=u.nm), None]
upper_corner = [SpectralCoord(280, unit=u.nm), None]
mg_crop = mg_ii.crop(lower_corner, upper_corner)

ax = mg_crop.plot()
plt.xlabel("Solar X")
plt.ylabel("Solar Y")

plt.show()

###############################################################################
# Imagine there's a really cool feature at (-338", 275"), how can you plot
# the spectrum at that location?

lower_corner = [None, SkyCoord(-338 * u.arcsec, 275 * u.arcsec, frame=Helioprojective)]
upper_corner = [None, SkyCoord(-338 * u.arcsec, 275 * u.arcsec, frame=Helioprojective)]
mg_ii.crop(lower_corner, upper_corner).plot()

plt.show()

###############################################################################
# Now, you may also be interested in knowing the time that was this observation taken.
# There is some information in `.meta`,
#
# But this is mostly about the observation in general.
# Times of individual scans are saved in ``.extra_coords['time']```.
# Getting access to them can be done in two ways:
#
# 1. Using axis_world_coords
# 2. "slice" the data object and access the .global_coords attribute.

# First let just see what we get from .meta
print(mg_ii.meta)

# Method 1: The first index is to escape the tuple that `axis_world_coords` returns
print(mg_ii.axis_world_coords("time", wcs=mg_ii.extra_coords)[0][0].isot)

# Method 2
print(mg_ii[0].global_coords["time"].isot)
