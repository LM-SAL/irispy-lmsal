from copy import copy

import astropy.modeling.models as m
import astropy.units as u
import gwcs
import gwcs.coordinate_frames as cf
import numpy as np
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS
from dkist.wcs.models import CoupledCompoundModel, VaryingCelestialTransform
from sunpy.coordinates.ephemeris import get_body_heliographic_stonyhurst
from sunpy.coordinates.frames import Helioprojective

from irispy.sji import SJICube
from irispy.utils import calculate_uncertainty
from irispy.utils.constants import BAD_PIXEL_VALUE_SCALED, BAD_PIXEL_VALUE_UNSCALED, DN_UNIT, READOUT_NOISE

__all__ = ["read_sji_lvl2"]


def _create_gwcs(hdulist: fits.HDUList) -> gwcs.WCS:
    """
    Creates the GWCS object for the SJI file.

    Parameters
    ----------
    hdulist : `astropy.io.fits.HDUList`
        The HDU list of the SJI file.

    Returns
    -------
    `gwcs.WCS`
        GWCS object for the SJI file.
    """
    pc_table = hdulist[1].data[:, hdulist[1].header["PC1_1IX"] : hdulist[1].header["PC2_2IX"] + 1].reshape(-1, 2, 2)
    crval_table = hdulist[1].data[:, hdulist[1].header["XCENIX"] : hdulist[1].header["YCENIX"] + 1]
    crpix = [hdulist[0].header["CRPIX1"], hdulist[0].header["CRPIX2"]]
    cdelt = [hdulist[0].header["CDELT1"], hdulist[0].header["CDELT2"]]
    celestial = VaryingCelestialTransform(
        crpix=crpix * u.pixel,
        cdelt=cdelt * u.arcsec / u.pixel,
        pc_table=pc_table * u.arcsec,
        crval_table=crval_table * u.arcsec,
    )
    base_time = Time(hdulist[0].header["STARTOBS"], format="isot", scale="utc")
    times = hdulist[1].data[:, hdulist[1].header["TIME"]] * u.s
    # We need to account for a non-zero time delta.
    base_time += times[0]
    times -= times[0]
    temporal = m.Tabular1D(
        np.arange(hdulist[1].data.shape[0]) * u.pix,
        lookup_table=times,
        fill_value=np.nan,
        bounds_error=False,
        method="linear",
    )
    forward_transform = CoupledCompoundModel("&", left=celestial, right=temporal)
    celestial_frame = cf.CelestialFrame(
        axes_order=(0, 1),
        unit=(u.arcsec, u.arcsec),
        reference_frame=Helioprojective(observer="earth", obstime=base_time),
        axis_physical_types=["custom:pos.helioprojective.lon", "custom:pos.helioprojective.lat"],
        axes_names=("Longitude", "Latitude"),
    )
    temporal_frame = cf.TemporalFrame(Time(base_time), unit=(u.s,), axes_order=(2,), axes_names=("Time (UTC)",))
    output_frame = cf.CompositeFrame([celestial_frame, temporal_frame])
    input_frame = cf.CoordinateFrame(
        axes_order=(0, 1, 2),
        naxes=3,
        axes_type=["PIXEL", "PIXEL", "PIXEL"],
        unit=(u.pix, u.pix, u.pix),
    )
    return gwcs.WCS(forward_transform, input_frame=input_frame, output_frame=output_frame)


def _create_wcs(hdulist):
    """
    This is required as occasionally we need a normal WCS instead of a gWCS due
    to compatibility issues.

    This has been set to have an Earth Observer at the time of the
    observation.
    """
    wcses = []
    base_time = Time(hdulist[0].header["STARTOBS"], format="isot", scale="utc")
    times = hdulist[1].data[:, hdulist[1].header["TIME"]] * u.s
    # We need to account for a non-zero time delta.
    base_time += times[0]
    times -= times[0]
    for i in range(hdulist[0].header["NAXIS3"]):
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
        header["DATE_OBS"] = (base_time + times[i]).isot
        location = get_body_heliographic_stonyhurst("Earth", header["DATE_OBS"])
        header["HGLN_OBS"] = location.lon.value
        header["HGLT_OBS"] = location.lat.value
        wcses.append(WCS(header))
    return wcses


def read_sji_lvl2(filename, *, uncertainty=False, memmap=False):
    """
    Reads a level 2 SJI FITS.

    Parameters
    ----------
    filename: `str`
        Filename to read.
    uncertainty : `bool`, optional
        If `True` (not the default), will compute the uncertainty for the data (slower and
        uses more memory). If ``memmap=True``, the uncertainty is never computed.
    memmap : `bool`, optional
        If `True` (not the default), will not load arrays into memory, and will only read from
        the file into memory when needed. This option is faster and uses a
        lot less memory. However, because FITS scaling is not done on-the-fly,
        the data units will be unscaled, not the usual data numbers (DN).

    Returns
    -------
    `irispy.sji.SJICube`
        The data cube, using a gWCS.
    """
    with fits.open(filename, memmap=memmap, do_not_scale_image_data=memmap) as hdulist:
        hdulist.verify("silentfix")
        extra_coords = [
            ("exposure time", 0, hdulist[1].data[:, hdulist[1].header["EXPTIMES"]] * u.s),
            ("obs_vrix", 0, hdulist[1].data[:, hdulist[1].header["OBS_VRIX"]] * u.m / u.s),
            ("ophaseix", 0, hdulist[1].data[:, hdulist[1].header["OPHASEIX"]] * u.arcsec),
            ("pztx", 0, hdulist[1].data[:, hdulist[1].header["PZTX"]] * u.arcsec),
            ("pzty", 0, hdulist[1].data[:, hdulist[1].header["PZTY"]] * u.arcsec),
            ("slit x position", 0, hdulist[1].data[:, hdulist[1].header["SLTPX1IX"]] * u.arcsec),
            ("slit y position", 0, hdulist[1].data[:, hdulist[1].header["SLTPX2IX"]] * u.arcsec),
            ("xcenix", 0, hdulist[1].data[:, hdulist[1].header["XCENIX"]] * u.arcsec),
            ("ycenix", 0, hdulist[1].data[:, hdulist[1].header["YCENIX"]] * u.arcsec),
        ]
        data = hdulist[0].data
        data_nan_masked = hdulist[0].data
        out_uncertainty = None
        if memmap:
            data_nan_masked[data == BAD_PIXEL_VALUE_UNSCALED] = 0
            mask = None
            scaled = False
            unit = DN_UNIT["SJI_UNSCALED"]
        elif not memmap:
            data_nan_masked[data == BAD_PIXEL_VALUE_SCALED] = np.nan
            mask = data_nan_masked == BAD_PIXEL_VALUE_SCALED
            scaled = True
            unit = DN_UNIT["SJI"]
            if uncertainty:
                out_uncertainty = calculate_uncertainty(data, READOUT_NOISE["SJI"], DN_UNIT["SJI"])
        else:
            raise ValueError(f"memmap={memmap} is not supported.")
        map_cube = SJICube(
            data_nan_masked,
            _create_gwcs(hdulist),
            uncertainty=out_uncertainty,
            unit=unit,
            meta=hdulist[0].header,
            mask=mask,
            scaled=scaled,
            _basic_wcs=_create_wcs(hdulist),
        )
        [map_cube.extra_coords.add(*extra_coord) for extra_coord in extra_coords]
    return map_cube
