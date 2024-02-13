import tarfile
from copy import copy
from pathlib import Path

import astropy.modeling.models as m
import astropy.units as u
import gwcs
import gwcs.coordinate_frames as cf
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.time import Time, TimeDelta
from astropy.wcs import WCS
from dkist.wcs.models import CoupledCompoundModel, VaryingCelestialTransform
from sunpy import log
from sunpy.coordinates import Helioprojective
from sunpy.coordinates.ephemeris import get_body_heliographic_stonyhurst

from irispy.spectrograph import Collection, SGMeta, SpectrogramCube, SpectrogramCubeSequence
from irispy.utils import calculate_uncertainty, tar_extract_all
from irispy.utils.constants import DN_UNIT, READOUT_NOISE, SLIT_WIDTH

__all__ = ["read_spectrograph_lvl2"]


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
        pc_table=pc_table * u.pixel,
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
        axis_physical_types=[
            "custom:pos.helioprojective.lon",
            "custom:pos.helioprojective.lat",
        ],
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


def _pc_matrix(lam, angle_1, angle_2):
    return angle_1, -1 * lam * angle_2, 1 / lam * angle_2, angle_1


def read_spectrograph_lvl2(
    filenames,
    *,
    spectral_windows=None,
    uncertainty=False,
    memmap=False,
    revert_v34=False,
):
    """
    Reads IRIS level 2 spectrograph FITS from an OBS into an
    `.IRISSpectrograph` instance.

    .. warning::

        This function does not handle tar files.
        Either extract them manually or use `irispy.io.read_files`.

    Parameters
    ----------
    filenames: `list` of `str` or `str`
        Filename of filenames to be read. They must all be associated with the same
        OBS number.
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
    revert_v34 : `bool`, optional.
        Will undo the data and WCS flipping made to V34 observations.
        Defaults to `False`.

    Returns
    -------
    `ndcube.NDCollection`
    """
    with fits.open(filenames[0], memmap=memmap, do_not_scale_image_data=memmap) as hdulist:
        v34 = bool(hdulist[0].header["OBSID"].startswith("34"))
        hdulist.verify("silentfix")
        windows_in_obs = np.array(
            [hdulist[0].header[f"TDESC{i}"] for i in range(1, hdulist[0].header["NWIN"] + 1)],
        )
        # If spectral_window is not set then get every window.
        # Else take the appropriate windows
        if not spectral_windows:
            spectral_windows_req = windows_in_obs
            window_fits_indices = range(1, len(hdulist) - 2)
        else:
            spectral_windows_req = [spectral_windows] if isinstance(spectral_windows, str) else spectral_windows
            spectral_windows_req = np.asarray(spectral_windows_req, dtype="U")
            window_is_in_obs = np.asarray([window in windows_in_obs for window in spectral_windows_req])
            if not all(window_is_in_obs):
                missing_windows = window_is_in_obs is False
                msg = f"Spectral windows {spectral_windows[missing_windows]} not in file {filenames[0]}"
                raise ValueError(msg)
            window_fits_indices = np.nonzero(np.in1d(windows_in_obs, spectral_windows))[0] + 1
        data_dict = {window_name: [] for window_name in spectral_windows_req}

    for filename in filenames:
        with fits.open(filename, memmap=memmap, do_not_scale_image_data=memmap) as hdulist:
            hdulist.verify("silentfix")
            # Extract axis-aligned metadata.
            times = Time(hdulist[0].header["STARTOBS"]) + TimeDelta(
                hdulist[-2].data[:, hdulist[-2].header["TIME"]],
                format="sec",
            )
            fov_center = SkyCoord(
                Tx=hdulist[-2].data[:, hdulist[-2].header["XCENIX"]],
                Ty=hdulist[-2].data[:, hdulist[-2].header["YCENIX"]],
                unit=u.arcsec,
                frame=Helioprojective,
            )
            obs_vrix = hdulist[-2].data[:, hdulist[-2].header["OBS_VRIX"]] * u.m / u.s
            ophaseix = hdulist[-2].data[:, hdulist[-2].header["OPHASEIX"]]
            exposure_times_fuv = hdulist[-2].data[:, hdulist[-2].header["EXPTIMEF"]] * u.s
            exposure_times_nuv = hdulist[-2].data[:, hdulist[-2].header["EXPTIMEN"]] * u.s
            for i, window_name in enumerate(spectral_windows_req):
                meta = SGMeta(
                    hdulist[0].header,
                    window_name,
                    data_shape=hdulist[window_fits_indices[i]].data.shape,
                )
                exposure_times = exposure_times_nuv
                DN_unit = DN_UNIT["NUV"]
                readout_noise = READOUT_NOISE["NUV"]
                if "FUV" in meta.detector:
                    exposure_times = exposure_times_fuv
                    DN_unit = DN_UNIT["FUV"]
                    readout_noise = READOUT_NOISE["FUV"]
                meta.add("exposure time", exposure_times, None, 0)
                meta.add("exposure FOV center", fov_center, None, 0)
                meta.add("observer radial velocity", obs_vrix, None, 0)
                meta.add("orbital phase", ophaseix, None, 0)
                # Sit-and-stare have a CDELT of 0 which causes issues in WCS.
                # In this case, set CDELT to a small number.
                header = copy(hdulist[window_fits_indices[i]].header)
                # Account for a slit offset (POFFYNUV (45) or POFFYFUV (34))
                idx = 34 if meta.spectral_band == "FUV" else 45
                header["CRVAL3"] -= hdulist[-2].data[:, idx].mean() * (SLIT_WIDTH.value / 2)
                if header["CDELT3"] == 0:
                    header["CDELT3"] = 1e-10
                    ang1, ang2, ang3, ang4 = _pc_matrix(
                        header["CDELT3"] / header["CDELT2"],
                        hdulist[-2].data[:, 20].mean(),
                        hdulist[-2].data[:, 22].mean(),
                    )
                    header["PC2_2"] = ang1
                    header["PC2_3"] = ang2
                    header["PC3_2"] = ang3
                    header["PC3_3"] = ang4
                try:
                    wcs = WCS(header)
                except Exception as e:  # NOQA: BLE001
                    msg = (
                        f"WCS failed to load while reading one step of the raster due to {e}"
                        "The loading will continue but this will be missing in the final cube. "
                        f"Spectral window: {window_name}, step {i} in file: {filename}"
                    )
                    log.warning(msg)
                    continue
                out_uncertainty = None
                data_mask = None
                if not memmap:
                    data_mask = hdulist[window_fits_indices[i]].data == -200.0
                if uncertainty:
                    out_uncertainty = calculate_uncertainty(
                        hdulist[window_fits_indices[i]].data,
                        readout_noise,
                        DN_UNIT,
                    )
                if v34 and not revert_v34:
                    data = np.flip(hdulist[window_fits_indices[i]].data, axis=0)
                    header["PC2_3"] = -header["PC2_3"]
                    header["PC3_2"] = -header["PC3_2"]
                    header["CDELT3"] = np.abs(header["CDELT3"])
                    header["CRPIX3"] = header["NAXIS3"] - header["CRPIX3"] + 1
                    wcs = WCS(header)
                else:
                    data = hdulist[window_fits_indices[i]].data
                cube = SpectrogramCube(
                    data,
                    wcs=wcs,
                    uncertainty=out_uncertainty,
                    unit=DN_unit,
                    meta=meta,
                    mask=data_mask,
                )
                cube.extra_coords.add("time", 0, times, physical_types="time")
                data_dict[window_name].append(cube)
    window_data_pairs = [
        (window_name, SpectrogramCubeSequence(data_dict[window_name], common_axis=0))
        for window_name in spectral_windows_req
    ]
    return Collection(window_data_pairs, aligned_axes=(0, 1, 2))
