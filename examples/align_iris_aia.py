"""
============================
Aligning IRIS SJI to SDO/AIA
============================

In this example we will show how to align a rolled IRIS dataset to SDO/AIA.

You can get IRIS data with co-aligned SDO data (and more) from https://iris.lmsal.com/search/
"""
# sphinx_gallery_thumbnail_number = 5

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import sunpy.map
from aiapy.calibrate import update_pointing
from astropy.coordinates import SkyCoord
from astropy.time import Time, TimeDelta
from sunpy.coordinates.frames import Helioprojective
from sunpy.net import Fido
from sunpy.net import attrs as a

from irispy.io import read_files
from irispy.obsid import ObsID
from irispy.utils.utils import _download_data

###############################################################################
# We start with getting the data.
# This is done by downloading the data from the IRIS archive.
#
# In this case, we will use requests as to keep this example self contained
# but using your browser will also work.

url = [
    "http://www.lmsal.com/solarsoft/irisa/data/level2_compressed/2014/09/19/20140919_051712_3860608353/iris_l2_20140919_051712_3860608353_SJI_2832_t000.fits.gz",
]
_download_data(url)
sji_filename = "iris_l2_20140919_051712_3860608353_SJI_2832_t000.fits.gz"

###############################################################################
# We will now open the slit-jaw imager (SJI) file we just downloaded.

sji_2832 = read_files(sji_filename)
# Printing will give us an overview of the file.
print(sji_2832)
# ``.meta`` contains the entire FITS header from the primary HDU.
# Since it is very long, we won't actually print it here.
# print(sji_2832.meta)

###############################################################################
# Can't remember what is OBSID 3860608353?
# **irispy-lmsal** has an utility function that will provide more information.

ObsID(sji_2832.meta["OBSID"])

###############################################################################
# Now, we will plot the SJI. By default, **irispy-lmsal** will a coordinate grid,
# and will color the spatial axes separately.You can disable by passing in
# ``bypass_formatting=True`` to the ``plot`` method.

sji_2832.plot()

plt.show()

###############################################################################
# We also have the option of going directly to an individual scan.

sji_cut = sji_2832[45]
print(sji_cut)

###############################################################################
# This dataset has a peculiarity: the observation has a 45 degree roll.
# The image does not have a 45 degree rotation because plotting shows the data
# in the way they are written in the file.
# We will a coordinate grid to make this clear.
# You can also change the axis labels and ticks if you so desire.
# `WCSAxes provides us an API we can use. <https://docs.astropy.org/en/stable/visualization/wcsaxes/index.html>`__

plt.figure()
# **irispy** will format the axes that if you want to avoid, you can disable:
ax = sji_cut.plot(bypass_formatting=True)
ax.set_xlabel("Helioprojective Longitude (Solar-X) [arcsec]")
ax.set_ylabel("Helioprojective Latitude (Solar-Y) [arcsec]")
ax.set_title(f"IRIS SJI {sji_2832.meta['TWAVE1']}")

ax.coords.grid(color="orange", linestyle="-")

lon = ax.coords[0]
lat = ax.coords[1]

lon.set_ticklabel_position("all")
lat.set_ticklabel_position("all")
lon.set_axislabel(ax.get_xlabel(), color="black")
lon.set_ticklabel("black")
lat.set_axislabel(ax.get_ylabel(), color="red")
lat.set_ticklabel("red")

plt.show()

###############################################################################
# For example, let us cut out the top sunspot.
# We need to specify the corners for the cut (note the coordinate order is
# the same as the plotted image). Be aware that crop works in the default N/S frame,
# so it will crop along those axes where as the data is rotated.
# You will also need to create a proper bounding box, with 4 corners.
#
# ``crop`` will return you the smallest bounding box which contains those 4 points
# which we can see when we overlay the points we give it.

# We need to create a coordinate frame for the IRIS data.
# While this is stored in the WCS, getting a coordinate frame is a little more involved.
sji_frame = Helioprojective(observer="earth", obstime=sji_cut.global_coords["Time (UTC)"])
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
# We will want to align the data to the AIA.
# First we will want to pick a timestamp during the observation.
#
# Lets us now find the SJI observation where the time is closest to 06:00 on 2014-09-19.

time_target = Time("2014-09-19T06:00:00")
time_index = np.abs(sji_2832.time - time_target).argmin()
time_stamp = sji_2832.time[time_index].isot
print(time_index, time_stamp)

###############################################################################
# The fact that it is rolled 45 degrees makes manual alignment tricky
# and will illustrate the power of working with WCS.
# We will download an AIA 170 nm image from the VSO.
# Once we have acquired it, we will need to use **aiapy** to prep this image.

search_results = Fido.search(
    a.Time(time_stamp, Time(time_stamp) + TimeDelta(1 * u.minute), near=time_stamp),
    a.Instrument.aia,
    a.Wavelength(1700 * u.AA),
)
files = Fido.fetch(search_results)

aia_map = sunpy.map.Map(files[0])
aia_map = update_pointing(aia_map)

# You don't need to register AIA images unless you need them aligned to other AIA images.
# otherwise you are degrading the data as the affine transform is not perfect.
# But it is the last step to get a level 1.5 image.
# from aiapy.calibrate import register
#
# aia_map = register(aia_map)

###############################################################################
# Now we will crop the AIA image to cover the active region seen by IRIS.

aia_bottom_left = SkyCoord(-850 * u.arcsec, -50 * u.arcsec, frame=aia_map.coordinate_frame)
aia_top_right = SkyCoord(-650 * u.arcsec, 150 * u.arcsec, frame=aia_map.coordinate_frame)
aia_sub = aia_map.submap(aia_bottom_left, top_right=aia_top_right)
aia_sub.plot()

plt.show()

###############################################################################
# Now let's plot the IRIS field of view on the AIA image.
# This IRIS data has no observer coordinate information
# **irispy-lmsal** will set this to be at Earth.
# This will allow us to transform from IRIS to any another observer.
#
# Using :meth:`~.draw_quadrangle`, drawing regions is straightforward.

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
    transform=ax.get_transform(sji_frame),
)

plt.show()

###############################################################################
# We have a green square showing the region of the IRIS observation.
# To work with both IRIS and AIA data, it helps if the image axes are aligned,
# and for this we need to rotate one of them. We can either rotate SDO/AIA to the
# IRIS frame, or vice-versa.
#
# We will rotate the AIA data, using the inverse rotation of the IRIS frame

aia_rot = aia_sub.rotate(rmatrix=np.matrix(sji_frame.wcs.pc[:-1, :-1]).I)

# Crop the AIA FOV to match IRIS.
bl = aia_rot.wcs.world_to_pixel(sji_frame.pixel_to_world(*(0, 0) * u.pix))
tr = aia_rot.wcs.world_to_pixel(sji_frame.pixel_to_world(*(sji_cut.data.shape[0], sji_cut.data.shape[1]) * u.pix))

fig = plt.figure()
ax1 = fig.add_subplot(1, 2, 1, projection=sji_frame)
aia_rotate_crop = aia_rot.submap(bl * u.pix, top_right=tr * u.pix)
aia_rotate_crop.plot(axes=ax1, autoalign=True)

ax2 = fig.add_subplot(1, 2, 2, projection=sji_cut)
sji_cut.plot(axes=ax2, cmap="irissji2832", vmin=0, vmax=4500)
ax2.set_title(f"IRIS SJI {sji_2832.meta['TWAVE1']}Å")
ax2.grid(color="w", ls=":")
ax2.set_xlabel(" ")
ax2.set_ylabel(" ")

plt.show()
