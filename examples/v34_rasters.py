"""
=============================
Dealing with IRIS v34 rasters
=============================

In this example we will show irispy-lmsal deals with a v34 dataset by default and how to
undo that if you so desire.

These v34 are scans which raster from west to east instead of the default east to west.
"""

import matplotlib.pyplot as plt
import pooch

import astropy.units as u
from astropy.coordinates import SpectralCoord

from irispy.io import read_files

###############################################################################
# We start by downloading a raster v34 dataset from the IRIS archive.
#
# In this case, we will use ``pooch`` to keep this example self-contained
# but using your browser will also work.
#
# This is a very large sparse 64-step raster with a FOV of 63"x175" for
# a total of 8 complete raster scans.

raster_filename = pooch.retrieve(
    "https://www.lmsal.com/solarsoft/irisa/data/level2_compressed/2025/03/28/20250328_225628_3400109360/iris_l2_20250328_225628_3400109360_raster.tar.gz",
    known_hash="fe890425c61a8d36e08806df19957ca264e62e13d7cda7e8d0a8f896ddd73db1",
)

###############################################################################
# We will now open the raster file we just downloaded.

raster = read_files(raster_filename, memmap=False)
# We will also undo the v34 handling and read the data as it is in the file.
# This will be used later to compare the results.
raster_unflipped = read_files(raster_filename, memmap=False, revert_v34=True)

# Printing will give us an overview of the file.
print(raster)

###############################################################################
# We will now pull the 'Mg II k 2796' spectral line data from the raster.

mg_ii_k = raster["Mg II k 2796"]
mg_ii_k_unflipped = raster_unflipped["Mg II k 2796"]
print(mg_ii_k)

###############################################################################
# So by default, irispy-lmsal will read the v34 data, flipping the data so that it
# is in the same orientation as normal IRIS data and adjust the WCS accordingly.

fig = plt.figure()
mg_ii_k.plot(fig=fig, clip_interval=(0.1, 99.9) * u.percent)

###############################################################################
# To confirm this we will plot spectroheliogram for Mg II k core wavelength.
# We can use the ``crop`` method to get this information, this will
# require a `~.SpectralCoord` object from `astropy.coordinates`.

# None, means that the axis is not cropped
# Note that this has to be in axis order
lower_corner = [SpectralCoord(280, unit=u.nm), None]
upper_corner = [SpectralCoord(280, unit=u.nm), None]
mg_spec_crop = mg_ii_k[0].crop(lower_corner, upper_corner)

fig = plt.figure(figsize=(6, 12))
ax = fig.add_subplot(121, projection=mg_spec_crop.wcs)
ax.set_title("v34 flipped")
mg_spec_crop.plot(axes=ax, plot_axes=["x", "y"])

###############################################################################
# Now let us compare the results with the unflipped data.

mg_spec_unflipped_crop = mg_ii_k_unflipped[0].crop(lower_corner, upper_corner)

ax2 = fig.add_subplot(122, projection=mg_spec_unflipped_crop.wcs)
ax2.set_title("v34 unflipped")
mg_spec_unflipped_crop.plot(axes=ax2, plot_axes=["x", "y"])
fig.tight_layout()

plt.show()

###############################################################################
# As you can see, the v34 data is flipped in the y-axis and the WCS is
# adjusted accordingly.
