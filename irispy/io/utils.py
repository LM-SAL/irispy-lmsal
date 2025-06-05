import tarfile
from pathlib import Path

from astropy.io import fits

from sunpy import log as logger

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
    list
        If this is a mixed collection of file data types
    `irispy.sji.SJICube`
        If a singular SJI file was passed in
    `irispy.spectrogram.SpectrogramCube`
        If a singular raster file was passed in
    `NDCollection`
        If a singular set of IRIS AIA aligned files were passed in.
    """
    from irispy.io.sji import read_sji_lvl2  # , read_aia_cube
    from irispy.io.spectrograph import read_spectrograph_lvl2

    if isinstance(filename, str | Path):
        filename = [filename]
    filename = sorted(filename)
    to_add = []
    to_remove = []
    for file in filename:
        if tarfile.is_tarfile(file):
            path = Path(file.replace(".tar.gz", ""))
            path.mkdir(parents=True, exist_ok=True)
            with tarfile.open(file, "r") as tar:
                tar.extractall(path, filter="data")
                to_add.extend([path / file for file in tar.getnames()])
                to_remove.append(file)
    filename.extend(to_add)
    for remove in to_remove:
        filename.pop(filename.index(remove))
    returns = []
    for file in filename:
        instrume = fits.getval(file, "INSTRUME")
        logger.debug(f"Processing file: {file} with instrume: {instrume}")
        try:
            if instrume == "IRIS" or instrume.startswith("AIA"):
                returns.append(read_sji_lvl2(file, memmap=memmap, uncertainty=uncertainty, **kwargs))
            elif instrume == "SPEC":
                returns.append(
                    read_spectrograph_lvl2(
                        file, spectral_windows=spectral_windows, memmap=memmap, uncertainty=uncertainty, **kwargs
                    )
                )
            else:
                logger.warning(f"INSTRUME: {instrume} was not recognized and not loaded")
        except Exception as e:
            if allow_errors:
                logger.warning(f"File {file} failed to load with {e}")
                continue
            raise
    if len(returns) == 1:
        return returns[0]
    return returns
