"""
========================
Overview of irispy-lmsal
========================

In this example, we will show how to use irispy-lmsal to open, crop and plot data.
We will finish with an example of aligning to AIA data.

You can get all IRIS data with co-aligned SDO data (and more) from https://iris.lmsal.com/search/
"""
# sphinx_gallery_thumbnail_number = 5

from copy import deepcopy

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import sunpy.map
from aiapy.calibrate import update_pointing
from astropy.coordinates import SkyCoord, SpectralCoord
from astropy.time import Time, TimeDelta
from astropy.wcs.utils import wcs_to_celestial_frame
from sunpy.coordinates.frames import Helioprojective
from sunpy.net import Fido
from sunpy.net import attrs as a

from irispy.io import read_sji_lvl2, read_spectrograph_lvl2
from irispy.obsid import ObsID
from irispy.utils.utils import _download_data

###############################################################################
# We start with getting the data. This is done by downloading the data from the IRIS archive.
#
# In this case, we will use requests as to keep this example self contained
# but using your browser will also work.
#
# Using the urls:
# http://www.lmsal.com/solarsoft/irisa/data/level2_compressed/2018/01/02/20180102_153155_3610108077/iris_l2_20180102_153155_3610108077_raster.tar.gz
# and
# http://www.lmsal.com/solarsoft/irisa/data/level2_compressed/2014/09/19/20140919_051712_3860608353/iris_l2_20140919_051712_3860608353_SJI_2832_t000.fits.gz

urls = [
    "http://www.lmsal.com/solarsoft/irisa/data/level2_compressed/2018/01/02/20180102_153155_3610108077/iris_l2_20180102_153155_3610108077_raster.tar.gz",
    "http://www.lmsal.com/solarsoft/irisa/data/level2_compressed/2014/09/19/20140919_051712_3860608353/iris_l2_20140919_051712_3860608353_SJI_2832_t000.fits.gz",
]
_download_data(urls)
raster_filename = "iris_l2_20180102_153155_3610108077_raster_t000_r00000.fits"
sji_filename = "iris_l2_20140919_051712_3860608353_SJI_2832_t000.fits.gz"

###############################################################################
# Working with IRIS spectrograph files
#
# Note that when ``memmap=True``, the data values are read from the FITS file
# directly without the scaling to Float32, the data values are no longer in DN,
# but in scaled integer units that start at −2$^{16}$/2.

raster = read_spectrograph_lvl2(raster_filename, memmap=True, uncertainty=False)
# Provide an overview of the data
print(raster)
# Will give us all the keys that corresponds to wavelengths.
print(raster.keys())

###############################################################################
# We will now access one wavelength window and plot it.

mg_ii = raster["Mg II k 2796"][0]

print(mg_ii)

mg_ii.plot()

plt.show()

###############################################################################
# We can index to get the first index and plot it.

print(mg_ii[0])

mg_ii[0].plot()

plt.show()

###############################################################################
# Or we can get the spectrum.

mg_ii[120, 200].plot()

plt.show()

###############################################################################
#  The default plots take the units and labels from the FITS WCS information,
# and often do not come in the most useful units (e.g. wavelengths in metres).
# We can read the wavelengths of the Mg window by calling axis_world_coords for wl
# (wavelength), and redo the plot with a better scale:

(mg_wave,) = mg_ii.axis_world_coords("wl")
mg_wave.to("m")

fig, ax = plt.subplots()
# Here we index the data directly instead of other methods to get access to the data and plot it "raw"
plt.plot(mg_wave.to("nm"), mg_ii.data[120, 200])

plt.show()

###############################################################################
#  When we use the underlying data directly, we lose all the metadata and WCS information,
# so we would suggest not doing it normally but there will be times you will need to do this.
# If you are unfamiliar with WCS, the following links are quite useful:
# The base comes from astropy: https://docs.astropy.org/en/stable/wcs/index.html
# The plotting makes use of WCSAxes: https://docs.astropy.org/en/stable/visualization/wcsaxes/index.html
# Some of the higher-level utils are via NDCube, e.g. coordinate transformations: https://docs.sunpy.org/projects/ndcube/en/v2.0.0rc2/coordinates.html.
# We can also improve on the default spectrogram plot by adjusting some options:

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
# Now we will plot spectroheliogram for Mg II k core wavelength.
# We can use a new feature in ndcube 2.0 called crop to get this information
# that will require a SkyCoord and SpectralCoord object from astropy.
# In this case, we only need to use a SpectralCoord.

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
#  Now, you may also be interested in knowing the time that was this observation taken.
# There is some information in .meta

print(mg_ii.meta)

###############################################################################
# But this is mostly about the observation in general.
# Times of individual scans are saved in .extra_coords['time'].
# Getting access to them can be done in two ways:
#
# 1. Using axis_world_coords
# 2. "slice" the data object and access the .global_coords attribute.

# The first index is to escape the tuple that `axis_world_coords` returns
print(mg_ii.axis_world_coords("time", wcs=mg_ii.extra_coords)[0][0].isot)

print(mg_ii[0].global_coords["time"].isot)

###############################################################################
# Working with IRIS Slit-Jaw Image files
#
# We will now open the slit-jaw image file we downloaded at the start.
sji_2832 = read_sji_lvl2(sji_filename)

print(sji_2832)
print(sji_2832.meta)

###############################################################################
# Can't remember what is OBSID 3860608353? irispy-lmsal
# has an utility function that will print out some more information

ObsID(sji_2832.meta["OBSID"])

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

lon.set_ticks(spacing=20 * u.arcsec, color="white", exclude_overlapping=True)
lat.set_ticks(spacing=20 * u.arcsec, color="white", exclude_overlapping=True)

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


sji_frame = wcs_to_celestial_frame(sji_2832.wcs)

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

(time_sji,) = sji_2832.axis_world_coords("time", wcs=sji_2832.extra_coords)
time_target = Time("2014-09-19T06:00:00.0")
time_index = np.abs(time_sji - time_target).argmin()
time_stamp = time_sji[time_index].isot
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

sji_corrected_wcs = deepcopy(sji_2832.wcs)
sji_corrected_wcs.wcs.dateobs = sji_2832[45].global_coords["time"].isot
sji_corrected_wcs.wcs.aux.hgln_obs = aia_map.observer_coordinate.lon.to_value(u.deg)
sji_corrected_wcs.wcs.aux.hglt_obs = aia_map.observer_coordinate.lat.to_value(u.deg)
sji_corrected_wcs.wcs.aux.rsun_ref = aia_map.observer_coordinate.rsun.to_value(u.m)

# We have to re-create the coordinate frame.
sji_frame_corrected = wcs_to_celestial_frame(sji_corrected_wcs)

# Remove time to avoid doing this call in all cells below
sji_corrected_wcs_drop = sji_corrected_wcs.dropaxis(-1)

###############################################################################
# Using the draw_quadrangle method, drawing WCS regions becomes much simpler
# but this only works with sunpy 3.1.

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
    transform=ax.get_transform(sji_corrected_wcs_drop),
)

plt.show()

###############################################################################
#  So now we have the green square showing the region of the IRIS observations.
# To work with both IRIS and AIA data, it helps if the image axes are aligned,
# and for this we need to rotate one of them. We can either rotate AIA to the
# IRIS frame, or vice-versa.
#
# We will rotate the AIA data, using the inverse rotation of the IRIS frame

aia_rot = aia_sub.rotate(rmatrix=np.matrix(sji_corrected_wcs.wcs.pc[:-1, :-1]).I)

# Crop the AIA FOV to match IRIS.
bl = aia_rot.wcs.world_to_pixel(sji_corrected_wcs_drop.pixel_to_world(*(0, 0) * u.pix))
tr = aia_rot.wcs.world_to_pixel(
    sji_corrected_wcs_drop.pixel_to_world(*(sji_cut.data.shape[0], sji_cut.data.shape[1]) * u.pix)
)

fig = plt.figure()
ax1 = fig.add_subplot(1, 2, 1, projection=sji_corrected_wcs_drop)
aia_rotate_crop = aia_rot.submap(bl * u.pix, top_right=tr * u.pix)
aia_rotate_crop.plot(axes=ax1, autoalign=True)

ax2 = fig.add_subplot(1, 2, 2, projection=sji_cut)
sji_cut.plot(axes=ax2, cmap="irissji2832", vmin=0, vmax=4500)
ax2.set_title(f"IRIS SJI {sji_2832.meta['TWAVE1']}Å")
ax2.grid(color="w", ls=":")
ax2.set_xlabel(" ")
ax2.set_ylabel(" ")

plt.show()
