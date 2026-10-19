"""
================
Plot a quick SJI
================

In this example we will show how to plot an SJI dataset.
"""

import matplotlib.pyplot as plt

from irispy.data import sample
from irispy.io import read_files

###############################################################################
# We start with the SJI 2832 sample file, which ``irispy.data.sample`` downloads on first use.
# To use your own data, pass the path of a Level 2 SJI file from the
# `IRIS data archive <https://www.lmsal.com/hek/hcr>`__ to ``read_files`` instead.

sji_filename = sample.SJI_2832

###############################################################################
# We will now open the data using a helper function which is designed to read
# all files from a single observation.

sji_2832 = read_files(sji_filename)

###############################################################################
# Printing will give us an overview of the SJI dataset.

print(sji_2832)

###############################################################################
# We will now plot the IRIS SJI data. The axis labels and ticks can be changed through the
# `WCSAxes API <https://docs.astropy.org/en/stable/visualization/wcsaxes/index.html>`__.

# ``get_animation`` lets Sphinx Gallery render the sequence as an animation.
animation = sji_2832.plot().get_animation()
plt.title(f"IRIS SJI {sji_2832.meta['TWAVE1']:.0f} Å", pad=25)
# ``animation.save("sji_2832.mp4")`` saves it as a movie, which needs FFmpeg.

###############################################################################
# Next, we will output a frame of the SJI as a sunpy Map.

sji_map = sji_2832.to_maps(0)
print(sji_map)

###############################################################################
# We can now plot the SJI Map using sunpy Map's plotting capabilities.

fig = plt.figure(layout="constrained")
ax = fig.add_subplot(projection=sji_map.wcs)
sji_map.plot(axes=ax)

plt.show()
