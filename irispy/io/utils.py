import sys
import tarfile
import warnings
from collections.abc import Iterable
from pathlib import Path

from astropy.io import fits
from astropy.table import Table
from sunpy import log

from irispy.io.sji import read_sji_lvl2
from irispy.io.spectrograph import read_spectrograph_lvl2

__all__ = ["fits_info", "read_files", "tar_extract_all"]


def tar_extract_all(filename: Path | str) -> list[Path]:
    """
    Extract all files a given tarfile and return the files.

    It will extract the files in the same folder as the tarfile.

    Parameters
    ----------
    filename : pathlib.Path, str
        The tarfile to extract.

    Returns
    -------
    list of pathlib.Path
        The list of files extracted from the tarfile.
    """
    filename = Path(filename)
    final_location = filename.parent / filename.stem
    final_location.mkdir(parents=True, exist_ok=True)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        with tarfile.open(filename, "r") as tar:
            tar.extractall(final_location)
            return [final_location / file for file in tar.getnames()]


def fits_info(filename: str) -> None:
    """
    Prints information about the extension of a raster or SJI level 2 data
    file.

    Parameters
    ----------
    filename : str
        Filename to load.
    """

    def get_description(idx, idx_mod):
        if idx == 0 and idx_mod == 0:
            return "Primary Header (no data)"
        if idx not in [num_extensions - 2, num_extensions - 1]:
            text_modifier = f" ({header[f'TDET{idx+idx_mod}'][:3]})" if "SPEC" in header[f"TDET{idx+idx_mod}"] else ""
            return f"{header[f'TDESC{idx+idx_mod}'].replace('_', ' ')} ({header[f'TWMIN{idx+idx_mod}']:.0f} - {header[f'TWMAX{idx+idx_mod}']:.0f} AA{text_modifier})"
        return "Auxiliary data"

    results = [
        f"Filename: {Path(filename).absolute()}",
        f"Observation: {fits.getval(filename, 'OBS_DESC')}",
        f"OBS ID: {fits.getval(filename, 'OBSID')}",
    ]
    table = Table(
        names=("No.", "Name", "Ver", "Type", "Cards", "Dimensions", "Format", "Description"),
        dtype=("S4", "S10", "S4", "S10", "S4", "S20", "S10", "S100"),
    )
    with fits.open(filename) as hdulist:
        hdu_info = hdulist.info(output=False)
        header = hdulist[0].header
        num_extensions = len(hdulist)
        for idx in range(num_extensions):
            description = get_description(idx, 0 if "SPEC" in header["INSTRUME"] else 1)
            hdu_info[idx] = hdu_info[idx][:-1] + (description,)
            table.add_row(list(map(str, hdu_info[idx])))

    sys.stdout.write("\n".join(results) + "\n")
    table.pprint()
    sys.stdout.flush()


def read_files(
    filepath: list[Path | str] | Path | str,
    *,
    spectral_windows: Iterable[str] | str | None = None,
    uncertainty: bool = False,
    memmap: bool = True,
) -> dict[str, list]:
    """
    A wrapper function to read a raster or SJI level 2 data file.

    You can provide:

    - one SJI file
    - one raster file
    - list of raster file
    - list of SJI files
    - list of mixed raster and SJI files

    If any of them are tar files, they will be extracted.
    Then they will append to the list of files to read at the end.

    This does not glob, so you cannot use wildcards, so you need to provide the
    full path to the files.

    The order of the files will be the same as the input list, minus the tar files.

    Parameters
    ----------
    filepath : `list of `str` or `pathlib.Path`, `str`, `pathlib.Path`
        Filename(s) to load.
    spectral_windows: iterable of `str` or `str`
        Spectral windows to extract from files.
        Default is None which means it will extract all spectral windows.
        Only used for raster files.
    uncertainty : `bool`, optional
        If `True` (not the default), will compute the uncertainty for the data (slower and uses more memory).
        If you set ``memmap`` to be True, the uncertainty is never computed.
    memmap : `bool`, optional
        If `True` (the default), will not load arrays into memory, and will only read from
        the file into memory when needed. This option is faster and uses a
        lot less memory. However, because FITS scaling is not done on-the-fly,
        the data units will be unscaled, not the usual data numbers (DN).

    Returns
    -------
    dict
        The dictionary contains the following keys:
        "raster"
            All loaded `irispy.spectrogram.SpectrogramCube` or `irispy.spectrogram.Collection`
        "sji"
            All loaded `irispy.sji.SJICube` or `irispy.spectrogram.SpectrogramCube`.
    """
    if isinstance(filepath, Path | str):
        filepath = [Path(filepath)]
    data_cubes = {"raster": [], "sji": []}
    tar_files = []
    for file in filepath:
        if tarfile.is_tarfile(file):
            tar_files.append(file)
            log.debug(f"Extracting {file}")
            filepath.remove(file)
            filepath.extend(tar_extract_all(file))
    for file in filepath:
        instrument = fits.getval(file, "INSTRUME")
        if instrument == "SJI":
            data_cubes["sji"].append(read_sji_lvl2(file, memmap=memmap, uncertainty=uncertainty))
        elif instrument == "SPEC":
            data_cubes["raster"].append(
                read_spectrograph_lvl2(
                    file,
                    spectral_windows=spectral_windows,
                    memmap=memmap,
                    uncertainty=uncertainty,
                )
            )
        else:
            msg = f"Unsupported instrument: {instrument}"
            raise ValueError(msg)
    return data_cubes
