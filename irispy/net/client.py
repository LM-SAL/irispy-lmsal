import numpy as np
import requests

import astropy.units as u
from astropy.time import Time

from sunpy.net import attrs as a
from sunpy.net.attr import AttrAnd, and_
from sunpy.net.base_client import BaseClient, QueryResponseTable
from sunpy.util.exceptions import warn_user

from irispy.net import attrs as iris_attrs

__all__ = ["LMSALClient"]

# The Heliophysics Coverage Registry (HCR) at LMSAL, which the IRIS search page also queries.
HCR_URL = "https://www.lmsal.com/hek/hcr"
# The HCR returns 400 observations unless asked for more, and at most 5000 (SolarSoft's iris_time2hcr.pro).
HCR_LIMIT = 5000


class LMSALQueryResponseTable(QueryResponseTable):
    # The other columns are hidden.
    display_keys = (
        "Start Time",
        "End Time",
        "OBSID",
        "Target",
        "NOAA AR",
        "XCen",
        "YCen",
        "FOV X",
        "FOV Y",
        "SDO",
        "SOT",
        "Product",
        "Wavelength",
        "Filename",
        "Size",
    )
    size_column = "Size"


def _quantity(value, unit):
    # Most of the HCR's numbers are single precision: drop the digits that the conversion to
    # double precision added.
    return None if value is None else float(f"{value:.7g}") * unit


def _tarball_size(event, key):
    # The HCR gives sizes in whole MiB, and -1 or nothing when there is no tarball.
    return _quantity(event.get(f"{key}_compsize_mb"), u.Mibyte) if f"{key}_compurl" in event else None


def _row(event, product, group, url, size):
    """
    The table row of one file, with `None` for what the HCR does not give.

    ``group`` is the HCR's entry for the file's product, empty for the AIA and SOT
    tarballs.
    """
    raster = product == "raster"
    fov = (event.get("raster_fovx"), event.get("raster_fovy")) if raster else (group.get("fovx"), group.get("fovy"))
    lines = event.get("spectral_window_names") if raster else None
    return {
        "Start Time": Time(group.get("start_time") or event["startTime"]),
        "End Time": Time(group.get("stop_time") or event["stopTime"]),
        "OBSID": int(event["obsId"]),
        "Target": event["target"],
        "NOAA AR": int(event["noaaNum"]) if event.get("noaaNum") else None,
        "XCen": round(event["xCen"], 1) * u.arcsec,
        "YCen": round(event["yCen"], 1) * u.arcsec,
        "FOV X": None if fov[0] is None else round(fov[0], 1) * u.arcsec,
        "FOV Y": None if fov[1] is None else round(fov[1], 1) * u.arcsec,
        "SDO": "martinaia_compurl" in event,
        "SOT": "martinsot_compurl" in event or "martinsp_compurl" in event,
        "Product": product,
        # The HCR's spectral_window_names lists the lines in a raster's windows, not the
        # windows, so it does not give the wavelength range of the raster.
        "Wavelength": _quantity(group.get("wavelnth"), u.AA),
        "Filename": url.rsplit("/", 1)[-1],
        # The HCR gives sizes in whole MiB.
        "Size": (np.nan if size is None else size) * u.Mibyte,
        "Observation Start": Time(event["startTime"]),
        "Title": event.get("obsTitle"),
        "OBS Description": event.get("goal"),
        "Science Objectives": event.get("sciObjectives"),
        "Planners": event.get("planners"),
        "Roll Angle": _quantity(event.get("roll_angle"), u.deg),
        "Duration": _quantity(event.get("duration_hours"), u.h),
        "Repeats": event.get("iris_repeats"),
        "Total Frames": event.get("iris_totframes"),
        "Images": group.get("num_images"),
        "Rasters": group.get("numrasters"),
        "Steps": group.get("raster_numsteps"),
        "Step Size": _quantity(group.get("raster_stepsize"), u.arcsec),
        "Exposure Time": _quantity(group.get("mean_exptime"), u.s),
        "Min Exposure Time": _quantity(group.get("min_exptime"), u.s),
        "Max Exposure Time": _quantity(group.get("max_exptime"), u.s),
        "Cadence": _quantity(group.get("cadence_avg_asrun"), u.s),
        "Planned Cadence": _quantity(group.get("cadence_avg_planned"), u.s),
        "Step Cadence": _quantity(group.get("basic_cadence_avg_asrun"), u.s),
        "Spatial Summing": event.get("spat_sum"),
        "FUV Spectral Summing": event.get("fuv_spec_sum"),
        "NUV Spectral Summing": event.get("nuvsji_spec_sum"),
        "Lines": lines and lines.strip("'"),
        "Data Fraction": event.get("data_fraction_found"),
        "Bad Data": event.get("iris_baddata"),
        "Rotation Tracking": event.get("rotation_tracking"),
        "SDO Size": _tarball_size(event, "martinaia"),
        "SOTFG Size": _tarball_size(event, "martinsot"),
        "SOTSP Size": _tarball_size(event, "martinsp"),
        "Movie URL": group.get("url_mp4"),
        "Thumbnail URL": group.get("url_thumb"),
        "Context URL": group.get("url_context"),
        "Timeline URL": event.get("timeline_gif_url"),
        "URL": url,
    }


class LMSALClient(BaseClient):
    """
    Searches the IRIS archive at LMSAL for Level 2 data and downloads it.

    Search with `a.Time <sunpy.net.attrs.Time>`, which finds the observations that overlap
    the time range, and ``a.Instrument("IRIS")`` or ``a.Provider("LMSAL")``; the
    `irispy.net.attrs` (``a.iris``) ``ObsID``, ``Description``, ``Target``, ``XCen``, ``YCen``
    and ``Product`` refine the search. ``a.Instrument("IRIS")`` gives the results of the IRIS
    archive and of the Virtual Solar Observatory, which lists IRIS under the same provider, so
    the VSO answers ``a.Provider("LMSAL")`` too; with an ``a.iris`` attr only the IRIS archive
    is searched, as the VSO does not know them. ``result["lmsal"]`` holds the IRIS archive's
    files.
    The archive answers about six observations a second, so a search of a month takes a minute.

    Each row of the result is one file: the raster tarball and the slit-jaw files of each
    observation, or with ``a.iris.Product`` only the files of that product, such as the
    co-aligned AIA cubes. Full-disk mosaics are not found; they are listed on the
    `IRIS mosaic page <https://iris.lmsal.com/mosaic.html>`__.

    The table shows each file's ``Start Time`` and ``End Time``, the observation's ``OBSID``,
    ``Target``, ``NOAA AR`` and pointing (``XCen``, ``YCen``), the file's field of view
    (``FOV X``, ``FOV Y``), whether the observation has co-aligned AIA or Hinode/SOT data
    (``SDO``, ``SOT``), the ``Product``, the slit-jaw ``Wavelength``, the ``Filename`` and the
    ``Size``. Its other columns are hidden: ``result["lmsal"].keys()`` lists every column,
    ``result["lmsal"].show("Filename", "Cadence")`` shows the columns named and
    ``result["lmsal"].show()`` all of them. They are:

    - The ``Observation Start``, ``Title``, ``OBS Description`` (such as
      ``"Large dense 320-step raster"``), ``Science Objectives`` and ``Planners`` of the
      observation.
    - Its ``Roll Angle``, ``Duration``, ``Repeats`` (how many times the OBS ran) and
      ``Total Frames`` (spectrograph exposures); ``Spatial Summing``, ``FUV Spectral Summing``
      and ``NUV Spectral Summing`` (of the NUV spectra and the slit-jaw images);
      ``Data Fraction``, the fraction of its data that the archive has, ``Bad Data`` and
      ``Rotation Tracking``.
    - The ``Images`` in a slit-jaw file; the ``Rasters`` in the raster tarball, their
      ``Steps`` and ``Step Size``, and ``Lines``, the lines in the raster's spectral windows.
    - ``Exposure Time`` (the mean, of slit-jaw files only), ``Min Exposure Time`` and
      ``Max Exposure Time``; ``Cadence``, the time between images or between rasters, and
      ``Planned Cadence``; ``Step Cadence``, the time between raster steps.
    - ``SDO Size``, ``SOTFG Size`` and ``SOTSP Size``, the sizes of the observation's AIA and
      Hinode/SOT tarballs.
    - The ``URL`` of the file, the ``Movie URL``, ``Thumbnail URL`` and ``Context URL``
      (slit-jaw files only) of its quick-look images, and the ``Timeline URL``, an image of
      the IRIS timeline of the day.

    The table masks what the archive does not give, such as the ``Wavelength`` of a raster
    tarball, the field of view of the AIA and SOT tarballs or the ``Repeats`` of early
    observations, and leaves out a column that none of the files found has.

    `Fido.fetch <sunpy.net.fido_factory.UnifiedDownloaderFactory.fetch>` puts the files of each
    observation in a folder named after it, such as ``20211001_060925_3683602040``, which
    `irispy.io.read_files` reads in one go. :ref:`irispy-tutorial-acquiring-data` shows a
    search and a download.

    Examples
    --------
    >>> from sunpy.net import Fido, attrs as a
    >>> import irispy.net
    >>> result = Fido.search(
    ...     a.Time("2021-10-01 06:00", "2021-10-01 07:00"),
    ...     a.Provider("LMSAL"),
    ...     a.iris.Target("AR"),
    ... )  # doctest: +REMOTE_DATA
    >>> print(
    ...     result["lmsal"]["OBSID", "Product", "Filename", "Size"]
    ... )  # doctest: +REMOTE_DATA
      OBSID    Product                         Filename                          Size
                                                                                Mibyte
    ---------- ------- -------------------------------------------------------- ------
    3683602040     sji iris_l2_20211001_060925_3683602040_SJI_1330_t000.fits.gz    3.0
    3683602040     sji iris_l2_20211001_060925_3683602040_SJI_1400_t000.fits.gz    3.0
    3683602040     sji iris_l2_20211001_060925_3683602040_SJI_2796_t000.fits.gz    5.0
    3683602040     sji iris_l2_20211001_060925_3683602040_SJI_2832_t000.fits.gz    6.0
    3683602040  raster         iris_l2_20211001_060925_3683602040_raster.tar.gz   83.0
    >>> files = Fido.fetch(result["lmsal"])  # doctest: +SKIP
    """

    @property
    def info_url(self):
        return "https://iris.lmsal.com/search/"

    def search(self, *query):
        """
        Search the IRIS archive.

        Parameters
        ----------
        *query : `sunpy.net.attr.Attr`
            The attrs to search with, combined with AND.

        Returns
        -------
        `sunpy.net.base_client.QueryResponseTable`
            One row per file.
        """
        query = and_(*query)
        params = {
            "cmd": "search-events3",
            "outputformat": "json",
            "instrument": "IRIS",
            "hasData": "true",
            "limit": HCR_LIMIT,
        }
        products = {"raster", "sji"}
        text = None
        for attr in query.attrs if isinstance(query, AttrAnd) else [query]:
            match attr:
                case a.Instrument() if str(attr.value).upper() == "IRIS":
                    pass
                case a.Provider() if str(attr.value).upper() == "LMSAL":
                    pass
                case a.Time():
                    params["startTime"], params["stopTime"] = attr.start.isot, attr.end.isot
                case iris_attrs.ObsID():
                    params["obsid"] = attr.value
                case iris_attrs.Target():
                    params["target"] = attr.value
                case iris_attrs.XCen():
                    params["minxcen"], params["maxxcen"] = attr.min.to_value(u.arcsec), attr.max.to_value(u.arcsec)
                case iris_attrs.YCen():
                    params["minycen"], params["maxycen"] = attr.min.to_value(u.arcsec), attr.max.to_value(u.arcsec)
                case iris_attrs.Description():
                    # The HCR's own text search (obsdesc) looks only at the title.
                    text = attr.value.lower()
                case iris_attrs.Product():
                    products = {attr.value}
                case _:
                    msg = f"LMSALClient cannot search with {attr!r}"
                    raise ValueError(msg)
        # A search of three months takes about 200 seconds.
        response = requests.get(HCR_URL, params=params, timeout=600)
        response.raise_for_status()
        events = response.json()["Events"]
        if len(events) == HCR_LIMIT:
            warn_user(f"The IRIS archive returns at most {HCR_LIMIT} observations; search a shorter time range.")
        rows = []
        for event in events:
            # Full-disk mosaics have no observation name.
            if "iris_obsshort" not in event:
                continue
            description = " ".join(event.get(key) or "" for key in ("obsTitle", "goal", "sciObjectives"))
            if text and text not in description.lower():
                continue
            files = [
                (
                    "raster"
                    if group["group_name"] == "Raster"
                    else "deconvolved"
                    if group["group_name"].endswith("_deconvolved")
                    else "sji",
                    group,
                    group.get("comp_data_url"),
                    group.get("comp_data_size_mb"),
                )
                for group in event["groups"]
            ]
            files += [
                ("sdo", {}, event.get("martinaia_compurl"), event.get("martinaia_compsize_mb")),
                ("sot", {}, event.get("martinsot_compurl"), event.get("martinsot_compsize_mb")),
                ("sot", {}, event.get("martinsp_compurl"), event.get("martinsp_compsize_mb")),
            ]
            rows += [
                _row(event, product, group, url, size)
                for product, group, url, size in files
                if product in products and url
            ]
        # The HCR returns the observations in no fixed order; a file's name starts with the
        # start of its observation.
        rows.sort(key=lambda row: row["Filename"])
        # Leave out what the HCR does not give, so that the table masks it, and the columns
        # that no row has.
        names = [key for key in rows[0] if any(row[key] is not None for row in rows)] if rows else None
        rows = [{key: value for key, value in row.items() if value is not None} for row in rows]
        return LMSALQueryResponseTable(rows=rows, names=names, client=self)

    def fetch(self, query_results, *, path, downloader, **kwargs):  # NOQA: ARG002
        """
        Queue the files of ``query_results`` on ``downloader``, each in the folder of
        its observation.

        `Fido.fetch <sunpy.net.fido_factory.UnifiedDownloaderFactory.fetch>` calls this.
        """
        for row in query_results:
            # The archive keeps each observation in a folder named after it.
            observation = row["URL"].rsplit("/", 2)[-2]
            filename = str(path).format(file=f"{observation}/{row['Filename']}", **row.response_block_map)
            downloader.enqueue_file(row["URL"], filename=filename)

    @classmethod
    def _can_handle_query(cls, *query):
        optional = {
            a.Instrument,
            a.Provider,
            iris_attrs.ObsID,
            iris_attrs.Description,
            iris_attrs.Target,
            iris_attrs.XCen,
            iris_attrs.YCen,
            iris_attrs.Product,
        }
        # As sunpy's GOES XRSClient does, the client's own attrs only refine a search for the
        # instrument or the provider.
        return (
            cls.check_attr_types_in_query(query, {a.Time}, optional)
            and any(isinstance(q, a.Instrument | a.Provider) for q in query)
            and all(str(q.value).upper() == "IRIS" for q in query if isinstance(q, a.Instrument))
            and all(str(q.value).upper() == "LMSAL" for q in query if isinstance(q, a.Provider))
        )

    @classmethod
    def _attrs_module(cls):
        return "iris", "irispy.net.attrs"

    @classmethod
    def register_values(cls):
        return {
            a.Instrument: [("IRIS", "Interface Region Imaging Spectrograph")],
            a.Provider: [("LMSAL", "Lockheed Martin Solar and Astrophysics Laboratory.")],
            iris_attrs.Product: list(iris_attrs._PRODUCTS.items()),
        }
