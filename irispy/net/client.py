import numpy as np
import requests

import astropy.units as u
from astropy.time import Time

from sunpy.net import attrs as a
from sunpy.net.attr import AttrAnd, and_
from sunpy.net.base_client import BaseClient, QueryResponseTable
from sunpy.util.exceptions import warn_user

from irispy.net import attrs as iris_attrs

__all__ = ["IRISClient"]

# The Heliophysics Coverage Registry (HCR) at LMSAL, which the IRIS search page also queries.
HCR_URL = "https://www.lmsal.com/hek/hcr"
# The HCR returns 400 observations unless asked for more, and at most 5000 (SolarSoft's iris_time2hcr.pro).
HCR_LIMIT = 5000


class IRISQueryResponseTable(QueryResponseTable):
    hide_keys = ("URL",)
    size_column = "Size"


class IRISClient(BaseClient):
    """
    Searches the IRIS archive at LMSAL for Level 2 data and downloads it.

    Search with `a.Time <sunpy.net.attrs.Time>`, which finds the observations that overlap
    the time range, and ``a.Instrument("IRIS")`` or the `irispy.net.attrs` (``a.iris``):
    ``ObsID``, ``Description``, ``Target``, ``XCen``, ``YCen`` and ``Product``.
    The archive answers about six observations a second, so a search of a month takes a minute.

    Each row of the result is one file: the raster tarball and the slit-jaw files of each
    observation, or with ``a.iris.Product`` only the files of that product, such as the
    co-aligned AIA cubes. The ``SDO`` and ``SOT`` columns say whether the observation has
    co-aligned AIA or Hinode/SOT data. Full-disk mosaics are not found; they are listed on the
    `IRIS mosaic page <https://iris.lmsal.com/mosaic.html>`__.

    `Fido.fetch <sunpy.net.fido_factory.UnifiedDownloaderFactory.fetch>` puts the files of each
    observation in a folder named after it, such as ``20211001_060925_3683602040``, which
    `irispy.io.read_files` reads in one go. :ref:`irispy-tutorial-acquiring-data` shows a
    search and a download.

    Examples
    --------
    An ``a.iris`` attr keeps the search to the IRIS archive; with only ``a.Instrument("IRIS")``
    the Virtual Solar Observatory is searched too.

    >>> from sunpy.net import Fido, attrs as a
    >>> import irispy.net
    >>> result = Fido.search(
    ...     a.Time("2021-10-01 06:00", "2021-10-01 07:00"), a.iris.Target("AR")
    ... )  # doctest: +REMOTE_DATA
    >>> print(
    ...     result["iris"]["OBSID", "Product", "Filename", "Size"]
    ... )  # doctest: +REMOTE_DATA
      OBSID    Product                         Filename                          Size
                                                                                Mibyte
    ---------- ------- -------------------------------------------------------- ------
    3683602040     sji iris_l2_20211001_060925_3683602040_SJI_1330_t000.fits.gz    3.0
    3683602040     sji iris_l2_20211001_060925_3683602040_SJI_1400_t000.fits.gz    3.0
    3683602040     sji iris_l2_20211001_060925_3683602040_SJI_2796_t000.fits.gz    5.0
    3683602040     sji iris_l2_20211001_060925_3683602040_SJI_2832_t000.fits.gz    6.0
    3683602040  raster         iris_l2_20211001_060925_3683602040_raster.tar.gz   83.0
    >>> files = Fido.fetch(result["iris"])  # doctest: +SKIP
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
                    msg = f"IRISClient cannot search with {attr!r}"
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
                    group.get("comp_data_url"),
                    group.get("comp_data_size_mb"),
                )
                for group in event["groups"]
            ]
            files += [
                ("sdo", event.get("martinaia_compurl"), event.get("martinaia_compsize_mb")),
                ("sot", event.get("martinsot_compurl"), event.get("martinsot_compsize_mb")),
                ("sot", event.get("martinsp_compurl"), event.get("martinsp_compsize_mb")),
            ]
            for product, url, size in files:
                if product in products and url:
                    rows.append(
                        {
                            "Start Time": Time(event["startTime"]),
                            "End Time": Time(event["stopTime"]),
                            "OBSID": int(event["obsId"]),
                            "Title": event.get("obsTitle") or "",
                            "Target": event["target"],
                            "XCen": round(event["xCen"], 1) * u.arcsec,
                            "YCen": round(event["yCen"], 1) * u.arcsec,
                            "SDO": "martinaia_compurl" in event,
                            "SOT": "martinsot_compurl" in event or "martinsp_compurl" in event,
                            "Product": product,
                            "Filename": url.rsplit("/", 1)[-1],
                            # The HCR gives sizes in whole MiB.
                            "Size": (np.nan if size is None else size) * u.Mibyte,
                            "URL": url,
                        }
                    )
        # The HCR returns the observations in no fixed order.
        rows.sort(key=lambda row: (row["Start Time"], row["Filename"]))
        return IRISQueryResponseTable(rows=rows, client=self)

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
            iris_attrs.ObsID,
            iris_attrs.Description,
            iris_attrs.Target,
            iris_attrs.XCen,
            iris_attrs.YCen,
            iris_attrs.Product,
        }
        # A search by time alone is left to the other clients.
        return (
            cls.check_attr_types_in_query(query, {a.Time}, optional)
            and any(not isinstance(q, a.Time) for q in query)
            and all(str(q.value).upper() == "IRIS" for q in query if isinstance(q, a.Instrument))
        )

    @classmethod
    def _attrs_module(cls):
        return "iris", "irispy.net.attrs"

    @classmethod
    def register_values(cls):
        return {
            a.Instrument: [("IRIS", "Interface Region Imaging Spectrograph")],
            iris_attrs.Product: list(iris_attrs._PRODUCTS.items()),
        }
