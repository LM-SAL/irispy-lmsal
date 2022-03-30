from copy import copy

import astropy.units as u
import numpy as np
from astropy.io import fits
from astropy.time import Time, TimeDelta
from astropy.wcs import WCS

from irispy.sji import IRISMapCube, IRISMapCubeSequence
from irispy.utils import calculate_uncertainty
from irispy.utils.constants import BAD_PIXEL_VALUE_SCALED, BAD_PIXEL_VALUE_UNSCALED, DN_UNIT, READOUT_NOISE

__all__ = ["read_sji_lvl2"]


def read_sji_lvl2(filename, uncertainty=False, memmap=False):
    """
    Read level 2 SJI FITS into an IRISMapCube instance.

    Parameters
    ----------
    filename: `str`
        Filename to read.
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
    `irispy.sji.IRISMapCubeSequence`
    """
    list_of_cubes = []
    with fits.open(filename, memmap=memmap, do_not_scale_image_data=memmap) as hdulist:
        hdulist.verify("silentfix")
        startobs = hdulist[0].header.get("STARTOBS")
        startobs = Time(startobs) if startobs else None
        endobs = hdulist[0].header.get("ENDOBS")
        endobs = Time(endobs) if endobs else None
        meta = {
            "TWAVE1": hdulist[0].header.get("TWAVE1"),
            "TELESCOP": hdulist[0].header.get("TELESCOP"),
            "INSTRUME": hdulist[0].header.get("INSTRUME"),
            "DATA_LEV": hdulist[0].header.get("DATA_LEV"),
            "OBSID": hdulist[0].header.get("OBSID"),
            "STARTOBS": startobs,
            "ENDOBS": endobs,
            "SAT_ROT": hdulist[0].header["SAT_ROT"] * u.deg,
            "NBFRAMES": hdulist[0].data.shape[0],
            "OBS_DESC": hdulist[0].header.get("OBS_DESC"),
            "FOVX": hdulist[0].header["FOVX"] * u.arcsec,
            "FOVY": hdulist[0].header["FOVY"] * u.arcsec,
            "XCEN": hdulist[0].header["XCEN"] * u.arcsec,
            "YCEN": hdulist[0].header["YCEN"] * u.arcsec,
        }
        times = Time(hdulist[0].header["STARTOBS"]) + TimeDelta(
            hdulist[1].data[:, hdulist[1].header["TIME"]], format="sec"
        )
        # The IRIS FITS files contain the per-exposure WCS metadata
        # (reference coordinate and PCij matrix) in an extension, while the primary
        # header has only values averaged over the observation
        for i in range(hdulist[0].header["NAXIS3"]):
            data = hdulist[0].data[i]
            data_nan_masked = hdulist[0].data[i]
            out_uncertainty = None
            if memmap:
                data_nan_masked[data == BAD_PIXEL_VALUE_UNSCALED] = 0
                mask = None
                scaled = False
                unit = DN_UNIT["SJI_UNSCALED"]
                out_uncertainty = None
            elif not memmap:
                data_nan_masked[data == BAD_PIXEL_VALUE_SCALED] = np.nan
                mask = data_nan_masked == BAD_PIXEL_VALUE_SCALED
                scaled = True
                unit = DN_UNIT["SJI"]
                if uncertainty:
                    out_uncertainty = calculate_uncertainty(data, READOUT_NOISE["SJI"], DN_UNIT["SJI"])
            else:
                raise ValueError(f"memmap={memmap} is not supported.")
            pztx = hdulist[1].data[i, hdulist[1].header["PZTX"]] * u.arcsec
            pzty = hdulist[1].data[i, hdulist[1].header["PZTY"]] * u.arcsec
            exposure_time = hdulist[1].data[i, hdulist[1].header["EXPTIMES"]] * u.s
            obs_vrix = hdulist[1].data[i, hdulist[1].header["OBS_VRIX"]] * u.m / u.s
            ophaseix = hdulist[1].data[i, hdulist[1].header["OPHASEIX"]] * u.arcsec
            slit_pos_x = hdulist[1].data[i, hdulist[1].header["SLTPX1IX"]] * u.arcsec
            slit_pos_y = hdulist[1].data[i, hdulist[1].header["SLTPX2IX"]] * u.arcsec
            global_coords = [
                ("time", "time", times[i]),
                ("pztx", "custom: PZTX", pztx),
                ("pzty", "custom: PZTY", pzty),
                ("obs_vrix", "custom: OBS_VRIX", obs_vrix),
                ("ophaseix", "custom: OPHASEIX", ophaseix),
                ("exposure time", "obs.exposure", exposure_time),
                ("slit x position", "custom: SLTPX1IX", slit_pos_x * u.pix),
                ("slit y position", "custom: SLTPX2IX", slit_pos_y * u.pix),
            ]
            header = copy(hdulist[0].header)
            header.pop("NAXIS3")
            header.pop("PC3_1")
            header.pop("PC3_2")
            header.pop("CTYPE3")
            header.pop("CUNIT3")
            header.pop("CRVAL3")
            header.pop("CRPIX3")
            header.pop("CDELT3")
            header["NAXIS"] = 2
            header["CRVAL1"] = hdulist[1].data[i, hdulist[1].header["XCENIX"]]
            header["CRVAL2"] = hdulist[1].data[i, hdulist[1].header["YCENIX"]]
            header["PC1_1"] = hdulist[1].data[0, hdulist[1].header["PC1_1IX"]]
            header["PC1_2"] = hdulist[1].data[0, hdulist[1].header["PC1_2IX"]]
            header["PC2_1"] = hdulist[1].data[0, hdulist[1].header["PC2_1IX"]]
            header["PC2_2"] = hdulist[1].data[0, hdulist[1].header["PC2_2IX"]]
            wcs = WCS(header)
            map_cube = IRISMapCube(
                data_nan_masked,
                wcs,
                uncertainty=out_uncertainty,
                unit=unit,
                meta=meta,
                mask=mask,
                scaled=scaled,
            )
            [map_cube.global_coords.add(*global_coord) for global_coord in global_coords]
            list_of_cubes.append(map_cube)
    return IRISMapCubeSequence(list_of_cubes, meta=meta, common_axis=None, times=times)
