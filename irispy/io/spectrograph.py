import logging
import tarfile
import textwrap
from copy import copy
from pathlib import Path

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.time import Time, TimeDelta
from astropy.wcs import WCS
from sunpy.coordinates import Helioprojective

from sunraster.meta import Meta, SlitSpectrographMetaABC

from irispy.spectrograph import IRISCollection, IRISSpectrogramCube, IRISSpectrogramCubeSequence
from irispy.utils import calculate_uncertainty
from irispy.utils.constants import DN_UNIT, READOUT_NOISE, SLIT_WIDTH, SPECTRAL_BAND

__all__ = ["read_spectrograph_lvl2"]


def _pc_matrix(lam, angle_1, angle_2):
    return angle_1, -1 * lam * angle_2, 1 / lam * angle_2, angle_1


def read_spectrograph_lvl2(filenames, spectral_windows=None, uncertainty=False, memmap=False):
    """
    Reads IRIS level 2 spectrograph FITS from an OBS into an
    `.IRISSpectrograph` instance.

    Parameters
    ----------
    filenames: `list` of `str` or `str`
        Filename of filenames to be read. They must all be associated with the same
        OBS number.
        If you provide a tar file, it will be extracted at the same location.
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

    Returns
    -------
    `ndcube.NDCollection`
    """
    if isinstance(filenames, list) and len(filenames) == 1:
        filenames = filenames[0]
    if isinstance(filenames, str):
        if tarfile.is_tarfile(filenames):
            parent = Path(filenames.replace(".tar.gz", "")).mkdir(parents=True, exist_ok=True)
            with tarfile.open(filenames, "r") as tar:
                tar.extractall(parent)
                filenames = [parent / file for file in tar.getnames()]
        else:
            filenames = [filenames]

    # Collecting the window observations
    with fits.open(filenames[0], memmap=memmap, do_not_scale_image_data=memmap) as hdulist:
        hdulist.verify("silentfix")
        windows_in_obs = np.array(
            [hdulist[0].header["TDESC{0}".format(i)] for i in range(1, hdulist[0].header["NWIN"] + 1)]
        )
        # If spectral_window is not set then get every window.
        # Else take the appropriate windows
        if not spectral_windows:
            spectral_windows_req = windows_in_obs
            window_fits_indices = range(1, len(hdulist) - 2)
        else:
            if isinstance(spectral_windows, str):
                spectral_windows_req = [spectral_windows]
            else:
                spectral_windows_req = spectral_windows
            spectral_windows_req = np.asarray(spectral_windows_req, dtype="U")
            window_is_in_obs = np.asarray([window in windows_in_obs for window in spectral_windows_req])
            if not all(window_is_in_obs):
                missing_windows = window_is_in_obs == False
                raise ValueError(
                    f"Spectral windows {spectral_windows[missing_windows]} not in file {filenames[0]}"
                )
            window_fits_indices = np.nonzero(np.in1d(windows_in_obs, spectral_windows))[0] + 1
        data_dict = dict([(window_name, list()) for window_name in spectral_windows_req])

    for filename in filenames:
        with fits.open(filename, memmap=memmap, do_not_scale_image_data=memmap) as hdulist:
            hdulist.verify("silentfix")
            # Extract axis-aligned metadata.
            times = Time(hdulist[0].header["STARTOBS"]) + TimeDelta(
                hdulist[-2].data[:, hdulist[-2].header["TIME"]], format="sec"
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
                meta = IRISSGMeta(
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
                idx = 45
                if meta.spectral_band == "FUV":
                    idx = 34
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
                except Exception as e:
                    logging.warning(
                        f"WCS failed to load while reading one step of the raster due to {e}"
                        " The loading will continue but this will be missing in the final cube."
                        f" Spectral window: {window_name}, step {i} in file: {filename}"
                    )
                    continue
                out_uncertainty = None
                data_mask = None
                if not memmap:
                    data_mask = hdulist[window_fits_indices[i]].data == -200.0
                if uncertainty:
                    out_uncertainty = calculate_uncertainty(
                        hdulist[window_fits_indices[i]].data, readout_noise, DN_UNIT
                    )
                cube = IRISSpectrogramCube(
                    hdulist[window_fits_indices[i]].data,
                    wcs=wcs,
                    uncertainty=out_uncertainty,
                    unit=DN_unit,
                    meta=meta,
                    mask=data_mask,
                )
                cube.extra_coords.add("time", 0, times, physical_types="time")
                data_dict[window_name].append(cube)
    window_data_pairs = [
        (window_name, IRISSpectrogramCubeSequence(data_dict[window_name], common_axis=0))
        for window_name in spectral_windows_req
    ]
    return IRISCollection(window_data_pairs, aligned_axes=(0, 1, 2))


class IRISSGMeta(Meta, metaclass=SlitSpectrographMetaABC):
    def __init__(self, header, spectral_window, **kwargs):
        super().__init__(header, **kwargs)
        spectral_windows = np.array([self["TDESC{0}".format(i)] for i in range(1, self["NWIN"] + 1)])
        window_mask = np.array([spectral_window in window for window in spectral_windows])
        if window_mask.sum() < 1:
            raise ValueError(
                "Spectral window not found. "
                f"Input spectral window: {spectral_window}; "
                f"Spectral windows in header: {spectral_windows}"
            )
        elif window_mask.sum() > 1:
            raise ValueError(
                "Spectral window must be unique. "
                f"Input spectral window: {spectral_window}; "
                f"Ambiguous spectral windows in header: {spectral_windows[window_mask]}"
            )
        self._iwin = np.arange(len(spectral_windows))[window_mask][0] + 1

    def __str__(self):
        return textwrap.dedent(
            f"""
                IRISMeta
                --------
                Observatory:     {self.observatory}
                Instrument:      {self.instrument}
                Detector:        {self.detector}
                Spectral Window: {self.spectral_window}
                Spectral Range:  {self.spectral_range}
                Spectral Band:   {self.spectral_band}
                Dimensions:      {self.dimensions}
                Date:            {self.date_reference}
                OBS ID:          {self.observing_mode_id}
                OBS Description: {self.observing_mode_description}
                """
        )

    def __repr__(self):
        return f"{object.__repr__(self)}\n{str(self)}"

    def _construct_time(self, key):
        val = self.get(key)
        if val is not None:
            val = Time(val, format="fits", scale="utc")
        return val

    @property
    def spectral_window(self):
        return self.get(f"TDESC{self._iwin}")

    @property
    def detector(self):
        return self.get(f"TDET{self._iwin}")

    @property
    def instrument(self):
        return self.get("INSTRUME")

    @property
    def observatory(self):
        return self.get("TELESCOP")

    @property
    def processing_level(self):
        return self.get("DATA_LEV")

    @property
    def distance_to_sun(self):
        return self.get("DSUN_OBS") * u.m

    @property
    def date_reference(self):
        return self._construct_time("DATE_OBS")

    @property
    def date_start(self):
        return self.date_reference

    @property
    def date_end(self):
        return self._construct_time("DATE_END")

    @property
    def observing_mode_id(self):
        return int(self.get("OBSID"))

    # ---------- IRIS-specific metadata properties ----------
    @property
    def dimensions(self):
        return self.shape.tolist()

    @property
    def observing_mode_description(self):
        return self.get("OBS_DESC")

    @property
    def observing_campaign_start(self):
        """
        Start time of observing campaign.
        """
        return self._construct_time("STARTOBS")

    @property
    def observing_campaign_end(self):
        """
        End time of observing mode.
        """
        return self._construct_time("ENDOBS")

    @property
    def observation_includes_SAA(self):
        """
        Whether IRIS passed through SAA during observations.
        """
        return bool(self.get("SAA"))

    @property
    def satellite_rotation(self):
        """
        Satellite roll from solar north.
        """
        return self.get("SAT_ROT") * u.deg

    @property
    def exposure_control_triggers_in_observation(self):
        """
        Number of times automatic exposure control triggered during observing
        campaign.
        """
        return self.get("AECNOBS")

    @property
    def exposure_control_triggers_in_raster(self):
        """
        Number of times automatic exposure control was triggered during this
        raster.
        """
        return self.get("AECNRAS")

    @property
    def number_raster_positions(self):
        """
        Number of positions in raster.
        """
        self.get("NRASTERP")

    @property
    def spectral_range(self):
        """
        The spectral range of the spectral window.
        """
        return [self.get(f"TWMIN{self._iwin}"), self.get(f"TWMAX{self._iwin}")] * u.AA

    @property
    def spectral_band(self):
        """
        The spectral band of the spectral window.
        """
        return SPECTRAL_BAND[self.spectral_window]

    @property
    def raster_fov_width_y(self):
        """
        Width of the field of view of the raster in the Y (slit) direction.
        """
        return self.get("FOVY") * u.arcsec

    @property
    def raster_fov_width_x(self):
        """
        Width of the field of view of the raster in the X (rastering)
        direction.
        """
        return self.get("FOVX") * u.arcsec

    @property
    def fov_center(self):
        """
        Location of the center of the field of view.
        """
        return SkyCoord(
            Tx=self.get("XCEN"),
            Ty=self.get("YCEN"),
            unit=u.arcsec,
            frame=Helioprojective,
        )

    @property
    def automatic_exposure_control_enabled(self):
        return bool(self.get("IAECFLAG"))

    @property
    def tracking_mode_enabled(self):
        return bool(self.get("TR_MODE"))

    @property
    def observatory_at_high_latitude(self):
        """
        Whether IRIS passed through high Earth latitude during observations.
        """
        return bool(self.get("HLZ"))

    @property
    def spatial_summing_factor(self):
        """
        Number of pixels summed together in the spatial (Y/slit) direction.
        """
        return self.get("SUMSPAT")

    @property
    def spectral_summing_factor(self):
        """
        Number of pixels summed together in the spectral direction.
        """
        if "fuv" in self.detector.lower():
            return self.get("SUMSPTRF")
        else:
            return self.get("SUMSPTRN")
