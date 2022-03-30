import tarfile
from pathlib import Path

from astropy.io import fits

__all__ = ["fitsinfo", "read_files"]


def fitsinfo(filename):
    """
    Prints information about the extension of a raster or SJI IRIS Level 2 data
    file.

    Parameters
    ----------
    filename : str
        Filename to load.
    """
    with fits.open(filename) as hdulist:
        hdulist.info()
        hdr = hdulist[0].header
        print("Observation description: ", hdr["OBS_DESC"], "\n")
        nwin = hdr["NWIN"]
        modifer = ""
        for i in range(nwin):
            print(f"Extension No. {i+1} stores data and header of {hdr[f'TDESC{i+1}']}: ", end="")
            if "SJI" not in hdr["TDET{}".format(i + 1)]:
                modifer = f" ({hdr[f'TDET{i+1}'][0:3]})"
            print(f"{hdr[f'TWMIN{i+1}']:.2f} - {hdr[f'TWMAX{i+1}']:.2f} AA" + modifer)


def read_files(filename, uncertainty=False, memmap=False):
    """
    A wrapper function to read a raster or SJI IRIS Level 2 data file.

    You can provide one SJI image or a one raster image or a list of raster images.

    If you mix raster and SJI images, the function will raise an error.

    Parameters
    ----------
    filename : `list of `str`, `str`
        Filename(s) to load.
        If given a string, will load that file.
        If given a list of strings, it will check they are all raster files and load them.
    uncertainty : `bool`, optional
        If `True` (not the default), will compute the uncertainty for the data (slower and
        uses more memory). If `memmap=True`, the uncertainty is never computed.
    memmap : `bool`, optional
        If `True` (not the default), will not load arrays into memory, and will only read from
        the file into memory when needed. This option is faster and uses a
        lot less memory. However, because FITS scaling is not done on-the-fly,
        the data units will be unscaled, not the usual data numbers (DN).

    Returns
    -------
    The corresponding `irispy.sji.IRISMapCube` or `irispy.spectrogram.IRISSpectrogramCube`.
    """
    from irispy.io.sji import read_sji_lvl2
    from irispy.io.spectrograph import read_spectrograph_lvl2

    if isinstance(filename, str):
        if tarfile.is_tarfile(filename):
            path = Path(filename.replace(".tar.gz", ""))
            path.mkdir(parents=True, exist_ok=True)
            with tarfile.open(filename, "r") as tar:
                tar.extractall(path)
                filename = [path / file for file in tar.getnames()]
        else:
            filename = [filename]

    intrume = fits.getval(filename[0], "INSTRUME")
    all_instrume = [fits.getval(f, "INSTRUME") for f in filename]
    if not all([intrume == i for i in all_instrume]):
        raise ValueError("You cannot mix raster and SJI files.")

    if intrume == "SJI":
        if len(filename) > 1:
            raise ValueError("You cannot load more than one SJI file at a time.")
        return read_sji_lvl2(filename[0], memmap=memmap, uncertainty=uncertainty)
    elif intrume == "SPEC":
        return read_spectrograph_lvl2(filename, memmap=memmap, uncertainty=uncertainty)
    else:
        raise ValueError(f"Unsupported instrument: {intrume}")
