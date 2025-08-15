"""
=======================
Radiometric Calibration
=======================

In this example we will show how to perform radiometric calibration on IRIS data.

IRIS data are given in counts or Data Number units (DN).
To convert these to a flux in physical units (e.g., :math:`erg s^{-1} sr^{-1} cm^{-2} Å^{-1}`)
one must perform a radiometric calibration.

Please refer to `ITN26 for more information on the calibration process <https://iris.lmsal.com/itn26/calibration.html>`__.
"""

from irispy.utils.response import get_latest_response

###############################################################################
# We will start by getting the latest IRIS response information.
# ``irispy`` only has support for the most up to date calibration.

iris_response = get_latest_response()

###############################################################################
# To convert the spectral units from DN to flux one must do the following calculation:
#
# .. math::
#
#    \mathrm{Flux}(\mathrm{erg}\: \mathrm{s}^{-1}\: \mathrm{cm}^{-2} \text{Å}^{-1}\: \mathrm{sr}^{-1}) = \mathrm{Flux}(\mathrm{DN}) \frac{E_\lambda \cdot \mathrm{DN2PHOT\_SG}}{A_\mathrm{eff} \cdot \mathrm{Pix}_{xy} \cdot \mathrm{Pix}_{\lambda} \cdot t_\mathrm{exp} \cdot W_\mathrm{slit}},
#
# where :math:`E_\lambda \equiv h \cdot c / \lambda` is the photon energy (in erg),
# :math:`DN2PHOT\_SG` is the number of photons per DN,
# :math:`A_\mathrm{eff}` is the effective area (in :math:`cm^{-2}`),
# :math:`Pix_{xy}` is the size of the spatial pixels in radians (e.g., multiply the spatial binning factor by :math:`\pi/(180\cdot3600\cdot6)`),
# :math:`Pix_{\lambda}` is the size of the spectral pixels in :math:`Å`,
# :math:`t_\mathrm{exp}` is the exposure time in seconds and
# :math:`W_\mathrm{slit}` is the slit width in radians (:math:`W_\mathrm{slit} \equiv \pi/(180\cdot3600\cdot3)`).
#
# This is a complex equation and requires careful attention to units.
# Within ``irispy``, there is a function called ``convert_dn_to_flux`` that handles this process.
