"""
====================
A quick and easy SJI
====================

In this example we will show how to plot a South Pole SJI dataset.

You can get IRIS data with co-aligned SDO data (and more) from https://iris.lmsal.com/search/
"""
import matplotlib.pyplot as plt
import pooch

from irispy.io import read_files

###############################################################################
# We start with getting the data.
# This is done by downloading the data from the IRIS archive.
#
# In this case, we will use ``pooch`` as to keep this example self contained
# but using your browser will also work.

sji_filename = pooch.retrieve(
    "https://www.lmsal.com/solarsoft/irisa/data/level2_compressed/2023/02/11/20230211_083601_3880012095/iris_l2_20230211_083601_3880012095_SJI_2832_t000.fits.gz",
    known_hash="a4b4f108ed67aded876bb167f24f093e8401d836f7deb66c97a36da8dd226064",
)

###############################################################################
# We will now open the slit-jaw imager (SJI) file we just downloaded.

sji_2832 = read_files(sji_filename)
# Printing will give us an overview of the file.
print(sji_2832)

###############################################################################
# We will now plot the IRIS SJI data.
# You can also change the axis labels and ticks if you so desire.
# `WCSAxes provides us an API we can use. <https://docs.astropy.org/en/stable/visualization/wcsaxes/index.html>`__

# Note that the .get_animation() is used to animate this example and is not required normally.
ax = sji_2832.plot().get_animation()
plt.xlabel("Helioprojective Longitude (Solar-X) [arcsec]")
plt.ylabel("Helioprojective Latitude (Solar-Y) [arcsec]")
plt.title(f"IRIS SJI {sji_2832.meta['TWAVE1']}", pad=20)

plt.show()
