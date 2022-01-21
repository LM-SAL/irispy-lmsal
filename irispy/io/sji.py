import astropy.units as u
import numpy as np
from astropy.io import fits
from astropy.time import Time, TimeDelta
from astropy.wcs import WCS

from irispy import utils
from irispy.sji import IRISMapCube, IRISMapCubeSequence
from irispy.utils.constants import BAD_PIXEL_VALUE_SCALED, BAD_PIXEL_VALUE_UNSCALED

__all__ = ["read_sji_lvl2"]


def read_sji_lvl2(filenames, uncertainty=True, memmap=False):
    """
    Read level 2 SJI FITS into an IRISMapCube instance.

    Parameters
    ----------
    filenames: `list` of `str` or `str`
        Filename or filenames to be read. They must all be associated with the same
        OBS number.
    uncertainty : `bool`, optional
        Default value is `True`.
        If `True`, will compute the uncertainty for the data (slower and
        uses more memory). If `memmap=True`, the uncertainty is never computed.
    memmap : `bool`, optional
        Default value is `False`.
        If `True`, will not load arrays into memory, and will only read from
        the file into memory when needed. This option is faster and uses a
        lot less memory. However, because FITS scaling is not done on-the-fly,
        the data units will be unscaled, not the usual data numbers (DN).

    Returns
    -------
    `irispy.sji.IRISMapCube` or `irispy.sji.IRISMapCubeSequence`
    """
    list_of_cubes = []
    if type(filenames) is str:
        filenames = [filenames]
    for filename in filenames:
        # Open a fits file
        hdulist = fits.open(filename, memmap=memmap, do_not_scale_image_data=memmap)
        hdulist.verify("fix")
        # Derive WCS, data and mask for NDCube from fits file.
        wcs = WCS(hdulist[0].header)
        data = hdulist[0].data
        data_nan_masked = hdulist[0].data
        if memmap:
            data_nan_masked[data == BAD_PIXEL_VALUE_UNSCALED] = 0
            mask = None
            scaled = False
            unit = utils.DN_UNIT["SJI_UNSCALED"]
            out_uncertainty = None
        elif not memmap:
            data_nan_masked[data == BAD_PIXEL_VALUE_SCALED] = np.nan
            mask = data_nan_masked == BAD_PIXEL_VALUE_SCALED
            scaled = True
            # Derive unit and readout noise from the detector
            unit = utils.DN_UNIT["SJI"]
            readout_noise = utils.READOUT_NOISE["SJI"]
            # Derive uncertainty of data for NDCube from fits file.
            if uncertainty:
                out_uncertainty = (
                    u.Quantity(
                        np.sqrt((data_nan_masked * unit).to(u.photon).value + readout_noise.to(u.photon).value ** 2),
                        unit=u.photon,
                    )
                    .to(unit)
                    .value
                )
            else:
                out_uncertainty = None
        # Derive exposure time from detector.
        exposure_times = hdulist[1].data[:, hdulist[1].header["EXPTIMES"]] * u.s
        # Derive extra coordinates for NDCube from fits file.
        times = Time(hdulist[0].header["STARTOBS"]) + TimeDelta(
            hdulist[1].data[:, hdulist[1].header["TIME"]], format="sec"
        )
        pztx = hdulist[1].data[:, hdulist[1].header["PZTX"]] * u.arcsec
        pzty = hdulist[1].data[:, hdulist[1].header["PZTY"]] * u.arcsec
        xcenix = hdulist[1].data[:, hdulist[1].header["XCENIX"]] * u.arcsec
        ycenix = hdulist[1].data[:, hdulist[1].header["YCENIX"]] * u.arcsec
        obs_vrix = hdulist[1].data[:, hdulist[1].header["OBS_VRIX"]] * u.m / u.s
        ophaseix = hdulist[1].data[:, hdulist[1].header["OPHASEIX"]] * u.arcsec
        slit_pos_x = hdulist[1].data[:, hdulist[1].header["SLTPX1IX"]] * u.arcsec
        slit_pos_y = hdulist[1].data[:, hdulist[1].header["SLTPX2IX"]] * u.arcsec
        extra_coords = [
            ("time", 0, times),
            ("pztx", 0, pztx),
            ("pzty", 0, pzty),
            ("xcenix", 0, xcenix),
            ("ycenix", 0, ycenix),
            ("obs_vrix", 0, obs_vrix),
            ("ophaseix", 0, ophaseix),
            ("exposure time", 0, exposure_times),
            ("slit x position", 0, slit_pos_x * u.pix),
            ("slit y position", 0, slit_pos_y * u.pix),
        ]
        # Extraction of meta for NDCube from fits file.
        startobs = hdulist[0].header.get("STARTOBS", None)
        startobs = Time(startobs) if startobs else None
        endobs = hdulist[0].header.get("ENDOBS", None)
        endobs = Time(endobs) if endobs else None
        meta = {
            "TELESCOP": hdulist[0].header.get("TELESCOP", None),
            "INSTRUME": hdulist[0].header.get("INSTRUME", None),
            "DATA_LEV": hdulist[0].header.get("DATA_LEV", None),
            "TWAVE1": hdulist[0].header.get("TWAVE1", None),
            "OBSID": hdulist[0].header.get("OBSID", None),
            "STARTOBS": startobs,
            "ENDOBS": endobs,
            "SAT_ROT": hdulist[0].header["SAT_ROT"] * u.deg,
            "NBFRAMES": hdulist[0].data.shape[0],
            "OBS_DESC": hdulist[0].header.get("OBS_DESC", None),
            "FOVX": hdulist[0].header["FOVX"] * u.arcsec,
            "FOVY": hdulist[0].header["FOVY"] * u.arcsec,
            "XCEN": hdulist[0].header["XCEN"] * u.arcsec,
            "YCEN": hdulist[0].header["YCEN"] * u.arcsec,
        }
        map_cube = IRISMapCube(
            data_nan_masked,
            wcs,
            uncertainty=out_uncertainty,
            unit=unit,
            meta=meta,
            mask=mask,
            scaled=scaled,
        )
        [map_cube.extra_coords.add(*extra_coord) for extra_coord in extra_coords]
        list_of_cubes.append(map_cube)
        hdulist.close()
    if len(filenames) == 1:
        return list_of_cubes[0]
    else:
        # In Sequence, all cubes must have the same Observation Identification.
        if np.any([cube.meta["OBSID"] != list_of_cubes[0].meta["OBSID"] for cube in list_of_cubes]):
            raise ValueError("Inputed files must have the same Observation Identification")
        # In Sequence, all cubes must have the same passband.
        if np.any([cube.meta["TWAVE1"] != list_of_cubes[0].meta["TWAVE1"] for cube in list_of_cubes]):
            raise ValueError("Inputed files must have the same passband")
        return IRISMapCubeSequence(list_of_cubes, meta=meta, common_axis=0)
