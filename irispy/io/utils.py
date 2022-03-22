from astropy.io import fits


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


def read_files(filename):
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

    Returns
    -------
    The corresponding `irispy.sji.IRISMapCube` or `irispy.spectrogram.IRISSpectrogramCube`.
    """
    from irispy.io.sji import read_sji_lvl2
    from irispy.io.sp import read_spectrograph_lvl2

    if isinstance(filename, str):
        filename = [filename]

    intrume = fits.getval(filename[0], "INSTRUME")
    all_instrume = [fits.getval(f, "INSTRUME") for f in filename]
    if not all([intrume == i for i in all_instrume]):
        raise ValueError("You cannot mix raster and SJI files.")

    if intrume == "SJI":
        return read_sji_lvl2(filename[0])
    elif intrume == "SPEC":
        return read_spectrograph_lvl2(filename)
    else:
        raise ValueError(f"Unsupported instrument: {intrume}")
