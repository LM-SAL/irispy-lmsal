"""
============================
Aligning IRIS SJI to SDO/AIA
============================
In this example we will show how to align a rolled IRIS dataset to SDO/AIA.
You can get IRIS data with co-aligned SDO data (and more) from https://iris.lmsal.com/search/
"""
# sphinx_gallery_thumbnail_number = 5

from copy import deepcopy

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import pooch
import sunpy.map
from aiapy.calibrate import update_pointing
from astropy.coordinates import SkyCoord
from astropy.time import Time, TimeDelta
from astropy.wcs.utils import wcs_to_celestial_frame
from sunpy.net import Fido
from sunpy.net import attrs as a

from irispy.io import read_files
from irispy.obsid import ObsID

###############################################################################
# We start with getting the data.
# This is done by downloading the data from the IRIS archive.
#
# In this case, we will use ``pooch`` as to keep this example self contained
# but using your browser will also work.

sji_filename = pooch.retrieve(
    "http://www.lmsal.com/solarsoft/irisa/data/level2_compressed/2014/09/19/20140919_051712_3860608353/iris_l2_20140919_051712_3860608353_SJI_2832_t000.fits.gz",
    known_hash=None,
)

###############################################################################
# We will now open the slit-jaw image file we downloaded.

sji_2832 = read_files(sji_filename)

print(sji_2832)
print(sji_2832.meta)

###############################################################################
# Can't remember what is OBSID 3860608353? irispy-lmsal
# has an utility function that will print out some more information

print(ObsID(sji_2832.meta["OBSID"]))

###############################################################################
# Now to plot it like the raster example.

sji_2832.plot()

plt.show()

###############################################################################
# We also have the option of going directly to a scan.

sji_cut = sji_2832[45]
print(sji_cut)

###############################################################################
#  We can make a more professional looking image like so

plt.figure()
ax = sji_cut.plot()
ax.grid(color="w", ls=":")
ax.set_xlabel("Helioprojective Longitude (Solar-X) [arcsec]")
ax.set_ylabel("Helioprojective Latitude (Solar-Y) [arcsec]")
ax.set_title(f"IRIS SJI {sji_2832.meta['TWAVE1']}")
plt.colorbar(ax=ax)

plt.show()

###############################################################################
#  This dataset has a peculiarity: the observation has a 45 degree roll.
# The image does not have a 45 degree rotation because plotting shows the data
# in the way they are written in the file.
# In the above, we added a coordinate grid to make this more clear.
# You can also change the axis labels and ticks if you so desire.
# `WCSAxes provides us an API you can find here <https://wcsaxes.readthedocs.io/en/latest/index.html>`__

plt.figure()
ax = sji_cut.plot()
ax.grid(color="w", ls=":")
ax.set_xlabel("Helioprojective Longitude (Solar-X) [arcsec]")
ax.set_ylabel("Helioprojective Latitude (Solar-Y) [arcsec]")
ax.set_title(f"IRIS SJI {sji_2832.meta['TWAVE1']}")
plt.colorbar(ax=ax)

lon = ax.coords[0]
lat = ax.coords[1]

lon.display_minor_ticks(True)
lat.display_minor_ticks(True)

lon.set_minor_frequency(20)
lat.set_minor_frequency(20)

lon.set_ticks(spacing=20 * u.arcsec, color="white")
lat.set_ticks(spacing=20 * u.arcsec, color="white")
lon.set_ticklabel(exclude_overlapping=True)
lat.set_ticklabel(exclude_overlapping=True)

plt.show()

###############################################################################
# Because the FITS files have the WCS coordinate information, we can use this
# to work in solar coordinates instead of pixels in the array.
# For example, let us cut out the top sunspot.
# We need to specify the corners for the cut (note the coordinate order is
# the same as the plotted image).Be aware that crop works in the default N/S frame,
# so it will crop along those axes where as the data is rotated.
# You will also need to create a proper bounding box, with 4 corners.
# Please be aware that crop is under heavy development in ndcube and
# might change its behavior in the future. It will also return you the smallest
# pixel box which contains those 4 points, which we can see when we overlay
# the points we give it. To start we will need the coordinate frame of the
# SJI observation to use for transforms and plotting later on.

sji_frame = wcs_to_celestial_frame(sji_2832[0].wcs)
bbox = [
    SkyCoord(-750 * u.arcsec, 90 * u.arcsec, frame=sji_frame),
    SkyCoord(-750 * u.arcsec, 95 * u.arcsec, frame=sji_frame),
    SkyCoord(-700 * u.arcsec, 90 * u.arcsec, frame=sji_frame),
    SkyCoord(-700 * u.arcsec, 95 * u.arcsec, frame=sji_frame),
]
sji_crop = sji_cut.crop(*bbox)
ax = sji_crop.plot(vmin=0, vmax=5000)
# Plot each corner of the box
[ax.plot_coord(coord, "o") for coord in bbox]
plt.xlabel("Solar X")
plt.ylabel("Solar Y")

plt.show()

###############################################################################
# Print the time of the current observation (``sji_cut`` object).

print(sji_cut.global_coords["time"])

###############################################################################
# Lets us now find the SJI observation where the time is closest to 06:00 on 2014-09-19.

time_target = Time("2014-09-19T06:00:00")
time_index = np.abs(sji_2832.time - time_target).argmin()
time_stamp = sji_2832.time[time_index].isot
print(time_index, time_stamp)

###############################################################################
# Aligning IRIS SJI with AIA
#
# Let us know assume we want to align this IRIS observation with an AIA image.
# The fact that it is rolled 45 degrees makes it even more interesting,
# difficult to do manually, and will illustrate the power of working with WCS.
# We will download an AIA 170 nm image from the VSO (or you can download it
# from the link at the top). Once we have acquired it, we will need to
# use aiapy + sunpy to prep and load this image.

search_results = Fido.search(
    a.Time(time_stamp, Time(time_stamp) + TimeDelta(1 * u.minute), near=time_stamp),
    a.Instrument.aia,
    a.Wavelength(1700 * u.AA),
)
files = Fido.fetch(search_results)

aia_map = sunpy.map.Map(files[0])
aia_map = update_pointing(aia_map)

# You don't need to register aia images unless you need them aligned to other aia images.
# otherwise you are degrading the data as the affine transform is not perfect.
# But it is the last step to get a 1.5 image.
# aia_map = register(aia_map)

###############################################################################
#  Now we will crop the AIA data to cover the active region.

aia_bottom_left = SkyCoord(-850 * u.arcsec, -50 * u.arcsec, frame=aia_map.coordinate_frame)
aia_top_right = SkyCoord(-650 * u.arcsec, 150 * u.arcsec, frame=aia_map.coordinate_frame)
aia_sub = aia_map.submap(aia_bottom_left, top_right=aia_top_right)
aia_sub.plot()

plt.show()


###############################################################################
# Now let's plot the IRIS field of view on the AIA image using the information
# from the WCS coordinates. This IRIS data has no observer coordinate information
# so we are just going to pretend it's at AIA for this example.
# This will allow us to transform from IRIS to SDO/AIA.

sji_corrected_wcs = deepcopy(sji_2832[0].wcs)
sji_corrected_wcs.wcs.dateobs = sji_2832.time[45].isot
sji_corrected_wcs.wcs.aux.hgln_obs = aia_map.observer_coordinate.lon.to_value(u.deg)
sji_corrected_wcs.wcs.aux.hglt_obs = aia_map.observer_coordinate.lat.to_value(u.deg)
sji_corrected_wcs.wcs.aux.rsun_ref = aia_map.observer_coordinate.rsun.to_value(u.m)

# We have to re-create the coordinate frame.
sji_frame_corrected = wcs_to_celestial_frame(sji_corrected_wcs)

###############################################################################
# Using the draw_quadrangle method, drawing WCS regions becomes much simpler.

fig = plt.figure()
ax = plt.subplot(projection=aia_sub)
aia_sub.plot()
aia_sub.draw_quadrangle(
    (0, 0) * u.pix,
    width=sji_cut.data.shape[1] * u.pix,
    height=sji_cut.data.shape[0] * u.pix,
    edgecolor="green",
    linestyle="--",
    linewidth=2,
    transform=ax.get_transform(sji_corrected_wcs),
)

plt.show()

###############################################################################
#  So now we have the green square showing the region of the IRIS observations.
# To work with both IRIS and AIA data, it helps if the image axes are aligned,
# and for this we need to rotate one of them. We can either rotate AIA to the
# IRIS frame, or vice-versa.
#
# We will rotate the AIA data, using the inverse rotation of the IRIS frame

aia_rot = aia_sub.rotate(rmatrix=np.array(np.matrix(sji_corrected_wcs.wcs.pc).I), missing=0)

# Crop the AIA FOV to match IRIS.
bl = aia_rot.wcs.world_to_pixel(sji_corrected_wcs.pixel_to_world(*(0, 0) * u.pix))
tr = aia_rot.wcs.world_to_pixel(
    sji_corrected_wcs.pixel_to_world(*(sji_cut.data.shape[0], sji_cut.data.shape[1]) * u.pix)
)

fig = plt.figure()
ax1 = fig.add_subplot(1, 2, 1, projection=sji_corrected_wcs)
aia_rotate_crop = aia_rot.submap(bl * u.pix, top_right=tr * u.pix)
aia_rotate_crop.plot(axes=ax1, autoalign=True)

ax2 = fig.add_subplot(1, 2, 2, projection=sji_cut)
sji_cut.plot(axes=ax2, cmap="irissji2832", vmin=0, vmax=4500)
ax2.set_title(f"IRIS SJI {sji_2832.meta['TWAVE1']}Å")
ax2.grid(color="w", ls=":")
ax2.set_xlabel(" ")
ax2.set_ylabel(" ")

plt.show()
