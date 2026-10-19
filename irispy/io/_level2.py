"""
Classify the files supplied with IRIS observations from their primary header.
"""

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
    Return the product of an IRIS level 2 file, or raise `ValueError`.

    Parameters
    ----------
    file : `str` or `pathlib.Path`
        The file, named in the error.
    header : `astropy.io.fits.Header` or `None`
        Its primary header; `None` for a file that is not FITS.
    product : `str`, optional
        ``"SJI"`` or ``"raster"`` when the file must be of that product.
    """
    found, level = _classify(header) if header is not None else (None, None)
    if found is None or level != 2 or product not in (None, found):
        expected = f"{product} file" if product else "file"
        msg = f"{file} is not an IRIS level 2 {expected}"
        raise ValueError(msg)
    return found
