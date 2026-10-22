"""
======================================
Fix the dust pixels in slit-jaw images
======================================

Dust on the slit-jaw detectors darkens the same detector pixels in every image, which
the IRIS pipeline sets to zero.
This example removes them and shows the difference.
"""

import matplotlib.pyplot as plt
import numpy as np

from irispy.data import sample
from irispy.io import read_files
from irispy.utils.dust import calculate_dust_mask, remove_dust

###############################################################################
# We start with the SJI 2832 sample file, which ``irispy.data.sample`` downloads on first use.
# To use your own data, pass the path of a Level 2 SJI file from the
# `IRIS data archive <https://www.lmsal.com/hek/hcr>`__ to ``read_files`` instead.

sji_2832 = read_files(sample.SJI_2832)

###############################################################################
# `~irispy.utils.dust.remove_dust` knows where the dust is for each channel.
# It replaces those pixels with the values of the neighboring images or,
# when the dust covers them too, with a local median.

clean_2832 = remove_dust(sji_2832)

###############################################################################
# To only mask the dust pixels and keep their values, use
# ``sji_2832.mask | calculate_dust_mask(sji_2832)`` as the cube's mask.

dust_mask = calculate_dust_mask(sji_2832)
print(f"{dust_mask[0].sum()} dust pixels in the first image")

###############################################################################
# Finally, we cut out a dust streak in the first image, compare it before and after,
# and show what changed.

original = sji_2832[0, 95:170, 350:430]
cleaned = clean_2832[0, 95:170, 350:430]
difference = cleaned - original.data * original.unit
limit = np.nanmax(np.abs(difference.data))
# The same stretch for both images, so that the dust pixels do not set it
vmin, vmax = np.nanpercentile(cleaned.data, [1, 99])

fig = plt.figure(figsize=(14, 5), layout="constrained")

ax0 = fig.add_subplot(131, projection=original.wcs)
original.plot(axes=ax0, vmin=vmin, vmax=vmax)
ax0.set_title("Original")

ax1 = fig.add_subplot(132, projection=cleaned.wcs)
cleaned.plot(axes=ax1, vmin=vmin, vmax=vmax)
ax1.set_title("Dust removed")

ax2 = fig.add_subplot(133, projection=difference.wcs)
difference.plot(axes=ax2, cmap="RdBu_r", vmin=-limit, vmax=limit)
ax2.set_title("Difference")
fig.colorbar(ax2.images[0], ax=ax2, label=str(difference.unit))

for ax in (ax1, ax2):
    ax.coords[1].set_ticklabel_visible(False)
    ax.coords[1].set_axislabel("")

plt.show()
