"""
Search attributes for `irispy.net.LMSALClient`, available as ``a.iris`` once
`irispy.net` is imported.

They narrow a search with ``a.Instrument("IRIS")`` or ``a.Provider("LMSAL")``.
"""

import astropy.units as u

from sunpy.net.attr import Range, SimpleAttr

__all__ = ["Description", "ObsID", "Product", "Target", "XCen", "YCen"]

_PRODUCTS = {
    "raster": "The raster tarball of the observation.",
    "sji": "The slit-jaw files of the observation.",
    "deconvolved": "The deconvolved slit-jaw files of the observation.",
    "sdo": "The tarball of AIA cubes co-aligned with the observation.",
    "sot": "The Hinode/SOT tarballs of the observation, SOTFG images and SOTSP maps.",
}


class ObsID(SimpleAttr):
    """
    The OBS ID of the observations, such as ``3683602040``.

    Parameters
    ----------
    value : `int` or `str`
    """

    def __init__(self, value):
        super().__init__(int(value))


class Description(SimpleAttr):
    """
    Text in the title, goal or science objectives of the observations, ignoring case.

    For example ``"Hinode"`` finds the observations made with Hinode and ``"sit-and-stare"`` the
    sit-and-stare observations.

    Parameters
    ----------
    value : `str`
    """


class Target(SimpleAttr):
    """
    The target of the observations as the archive labels it.

    For example ``"AR"`` (active region), ``"QS"`` (quiet Sun), ``"CH"`` (coronal hole) or
    ``"FI"`` (filament); the ``Target`` column of a search shows the labels.

    Parameters
    ----------
    value : `str`
    """

    def __init__(self, value):
        # The HCR matches the target case-sensitively: target=ar finds nothing.
        super().__init__(value.upper())


class _Pointing(Range):
    @u.quantity_input
    def __init__(self, min_: u.arcsec, max_: u.arcsec):
        super().__init__(min_, max_)

    def collides(self, other):
        return isinstance(other, self.__class__)


class XCen(_Pointing):
    """
    The solar x of the observations' pointing, from ``min_`` to ``max_``.

    Parameters
    ----------
    min_, max_ : `~astropy.units.Quantity`
        Angles from disk center, positive to the west.
    """


class YCen(_Pointing):
    """
    The solar y of the observations' pointing, from ``min_`` to ``max_``.

    Parameters
    ----------
    min_, max_ : `~astropy.units.Quantity`
        Angles from disk center, positive to the north.
    """


class Product(SimpleAttr):
    """
    The files to find in place of the raster tarball and slit-jaw files.

    One of ``"raster"``, ``"sji"``, ``"deconvolved"`` (deconvolved slit-jaw files), ``"sdo"``
    (co-aligned AIA cubes) or ``"sot"`` (Hinode/SOT tarballs); combine products with ``|``,
    such as ``Product("raster") | Product("sdo")``.
    `irispy.io.read_files` on a folder that holds both the slit-jaw and the deconvolved files
    gives two cubes for each band.

    Parameters
    ----------
    value : `str`
    """

    def __init__(self, value):
        super().__init__(value.lower())
