from pathlib import Path

from astropy.io import fits

from ndcube import NDCollection
from sunpy import log as logger

from irispy.io.sji import read_sji_lvl2
from irispy.io.spectrograph import read_spectrograph_lvl2

__all__ = ["fitsinfo", "read_files"]


def fitsinfo(filename):
    """
    Prints information about the extension of a raster or SJI level 2 data
    file.

    Parameters
    ----------
    filename : str
        Filename to load.
    """
    with fits.open(filename) as hdulist:
        hdulist.info()
        hdr = hdulist[0].header
        msg = f"Observation description: {hdr['OBS_DESC']}"
        logger.info(msg)
        modifier = ""
        for i in range(hdr["NWIN"]):
            msg = f"Extension No. {i + 1} stores data and header of {hdr[f'TDESC{i + 1}']}: "
            logger.info(msg)
            if "SJI" not in hdr[f"TDET{i + 1}"]:
                modifier = f" ({hdr[f'TDET{i + 1}'][:3]})"
            msg = f"{hdr[f'TWMIN{i + 1}']:.2f} - {hdr[f'TWMAX{i + 1}']:.2f} AA{modifier}"
            logger.info(msg)


def read_files(filename, *, spectral_windows=None, uncertainty=False, memmap=False, allow_errors=False, **kwargs):
    """
    A wrapper function to read any number of raster, SJI or IRIS-algined AIA
    data files.

    Parameters
    ----------
    filename : `list of `str`, `str`, `pathlib.Path`
        Filename(s) to load.
        If given a string, will load that file.
        If given a list of strings, load them.
    spectral_windows: iterable of `str` or `str`
        Spectral windows to extract from files. Default=None, implies, extract all
        spectral windows.
    uncertainty : `bool`, optional
        If `True` (not the default), will compute the uncertainty for the data (slower and
        uses more memory). If `memmap=True`, the uncertainty is never computed.
    memmap : `bool`, optional
        If `True` (not the default), will not load arrays into memory, and will only read from
        the file into memory when needed. This option is faster and uses a
        lot less memory. However, because FITS scaling is not done on-the-fly,
        the data units will be unscaled, not the usual data numbers (DN).
    allow_errors : `bool`, optional
        Will continue loading the files if one fails to load.
        Defaults to `False`.
    kwargs : `dict`, optional
        Additional keyword arguments to pass to the reader functions.

    Returns
    -------
    `NDCollection`
        With keys being the value of TDESC1, the values being the cube.
    """
    if isinstance(filename, (str, Path)):
        filename = [Path(filename)]
    filename = sorted(filename)
    filename = [Path(f) for f in filename]
    returns = {}
    for file in filename:
        instrume = ""
        describe = ""
        if file.name.endswith((".fits", ".fits.gz")):
            instrume = fits.getval(file, "INSTRUME")
            describe = fits.getval(file, "TDESC1")
        logger.debug(f"Processing file: {file} with instrume: {instrume}")
        try:
            if instrume in ["IRIS", "SJI"] or instrume.startswith("AIA"):
                returns[f"{describe}"] = read_sji_lvl2(file, memmap=memmap, uncertainty=uncertainty, **kwargs)
            # Only raster files are tar.gz files, so we lump them together
            # read_spectrograph_lvl2 itself handles tar.gz files.
            elif file.name.endswith(".tar.gz") or instrume == "SPEC":
                returns[f"{describe}"] = read_spectrograph_lvl2(
                    file, spectral_windows=spectral_windows, memmap=memmap, uncertainty=uncertainty, **kwargs
                )
            else:
                logger.warning(f"INSTRUME: {instrume} was not recognized and not loaded")
        except Exception as e:
            if allow_errors:
                logger.warning(f"File {file} failed to load with {e}")
                continue
            raise
    return NDCollection(returns.items()) if len(returns) > 1 else next(iter(returns.values()))
