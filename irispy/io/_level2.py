"""
Classify the files supplied with IRIS observations from their primary header.
"""

import os
from pathlib import Path

__all__ = []


def _classify(header):
    """
    Return ``(product, level)`` for the primary header of an IRIS data file.

    Parameters
    ----------
    header : `astropy.io.fits.Header`
        The primary header.

    Returns
    -------
    `tuple`
        ``product`` is ``"SJI"`` for slit-jaw images and the AIA and Hinode/SOT cubes
        supplied with IRIS observations, ``"raster"`` for spectrograph files and `None`
        for anything else, such as a full-disk mosaic. ``level`` is the data level from
        ``DATA_LEV`` (``LVL_NUM`` when missing), `None` when the header has neither.
    """
    instrume = str(header.get("INSTRUME", ""))
    if instrume == "SPEC":
        product = "raster"
    elif instrume in ("SJI", "IRIS") or instrume.startswith(("AIA", "SOT")):
        # Level 3 files also have INSTRUME="IRIS"; the level tells them apart.
        product = "SJI"
    else:
        product = None
    return product, header.get("DATA_LEV", header.get("LVL_NUM"))


def _check_level2(file, header, product=None):
    """
    Return the product of an IRIS level 2 file, or raise `ValueError`, or the `OSError`
    of a path that cannot be read.

    Parameters
    ----------
    file : `str`, `pathlib.Path` or `object`
        What was read, named in the error: a path as given, anything else by its type.
    header : `astropy.io.fits.Header` or `None`
        Its primary header; `None` for a file that is not FITS.
    product : `str`, optional
        ``"SJI"`` or ``"raster"`` when the file must be of that product.
    """
    if header is None and isinstance(file, (str, os.PathLike)) and not Path(file).is_dir():
        # A missing or unreadable path raises what opening it raises, as in irispy 0.9.1
        Path(file).open("rb").close()
    found, level = _classify(header) if header is not None else (None, None)
    if found is None or level != 2 or product not in (None, found):
        expected = f"{product} file" if product else "file"
        name = file if isinstance(file, (str, os.PathLike)) else f"The data given ({type(file).__name__})"
        msg = f"{name} is not an IRIS level 2 {expected}"
        raise ValueError(msg)
    return found
