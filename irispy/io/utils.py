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


def load(filename):
    """
    Temp workaround for a "unified" loader.

    Parameters
    ----------
    filename : str
        Filename to load.

    Returns
    -------
    NDCube
    """
    from sunraster.instr.iris import read_iris_spectrograph_level2_fits

    from .sji import read_sji_lvl2

    try:
        return read_sji_lvl2(filename)
    except Exception:
        try:
            return read_iris_spectrograph_level2_fits(filename)
        except Exception:
            print("Failed to load {filename}")
