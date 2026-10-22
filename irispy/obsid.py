import importlib.resources

import numpy as np

from astropy import units as u
from astropy.io import ascii as ascii_io

__all__ = ["ObsID"]


class ObsID(dict):
    """
    A class to convert the IRIS OBS ID to human-readable format.

    Parameters
    ----------
    obsid : `int` or `str`
        IRIS OBS ID to query, such as ``cube.meta.observing_mode_id`` of a cube or the ``OBSID`` header
        keyword. Needs to be a valid OBS ID, with 10 digits.

    Notes
    -----
    The keys and the types of their values:

    * ``"obsid"`` (`int`): the OBS ID.
    * ``"raster_desc"`` and ``"raster_fulldesc"`` (`str`): the description, such as
      ``"Large dense 96-step raster"``, and the same with the field of view and cadence.
    * ``"raster_fov"`` and ``"sji_fov"`` (`str`): the fields of view in arcsec, such as ``"31.35x120"``.
    * ``"raster_step"`` (`~astropy.units.Quantity`): the raster step, one per sub-raster for
      sparse/dense/sparse rasters.
    * ``"exptime"`` (`~astropy.units.Quantity`): the exposure time.
    * ``"spec_cadence"`` (`~astropy.units.Quantity`): the nominal cadence of the OBS table, which counts
      1 s per raster step; the cadence of the data is longer, as the cube's times show.
    * ``"linelist"``, ``"compression"``, ``"sji_cadence"``, ``"fuv_binning"``, ``"binning"`` and
      ``"sjis"`` (`str`): the options of the OBS table, such as ``"Flare linelist 1"``.
    * ``"readout"`` (`str`): ``"Non-simultaneous readout"`` or ``"Simultaneous readout"``, only in OBS
      table versions 34 and 36.

    ``options`` holds, for each option, the choices of the OBS table.

    OBS IDs from the OBS table versions 34, 36, 38 and 40 are supported.
    OBS IDs of 4200000000 and above are calibration and special sequences
    with no OBS table; the ``OBS_DESC`` and ``OBSTITLE`` header keywords
    describe them.

    Examples
    --------
    Quickly show OBS ID parameters

    >>> from irispy import obsid
    >>> obsid.ObsID(3675508564)
    IRIS OBS ID 3675508564
    ----------------------
    Description:            Large dense 96-step raster 31.35x120 96s
    SJI filters:                                    C II   Mg II w s
    SJI field of view:                                       120x120
    Exposure time:                                             8.0 s
    Binning:                               Spatial x 1, Spectral x 1
    FUV binning:                         FUV spectrally rebinned x 4
    SJI cadence:                                 SJI cadence default
    Compression:                                Lossless compression
    Linelist:                                       Flare linelist 1

    The data can be accessed as in a dictionary:

    >>> data = obsid.ObsID(3675508564)
    >>> data["exptime"]
    <Quantity 8. s>
    >>> data["linelist"]
    'Flare linelist 1'
    >>> data["raster_fov"]
    '31.35x120'
    >>> data["raster_step"]
    <Quantity 0.33 arcsec>
    >>> data["spec_cadence"]
    <Quantity 96. s>

    Sparse/dense/sparse rasters have one raster step per sub-raster:

    >>> obsid.ObsID(4191006128)["raster_step"]
    <Quantity [1.  , 0.33, 1.  ] arcsec>
    """

    def __init__(self, obsid) -> None:
        self.obsid = int(obsid)
        data, options = self._read_obsid(self.obsid)
        super().__init__(data)
        self.options = options

    def __repr__(self) -> str:
        return (
            "IRIS OBS ID {obsid}\n"
            "----------------------\n"
            "Description:       {raster_fulldesc:>45}\n"
            "SJI filters:       {sjis:>45}\n"
            "SJI field of view: {sji_fov:>45}\n"
            "Exposure time:     {exptime!s:>45}\n"
            "Binning:           {binning:>45}\n"
            "FUV binning:       {fuv_binning:>45}\n"
            "SJI cadence:       {sji_cadence:>45}\n"
            "Compression:       {compression:>45}\n"
            "Linelist:          {linelist:>45}"
        ).format(**self)

    @staticmethod
    def _exptime_to_quant(exptime):
        """
        Converts an 'exptime' string (used in IRIS tables and OBS_DESC) to a Quantity
        instance in seconds.
        """
        if exptime == "Exposure 1s":
            return 1.0 * u.s
        return float(exptime.split(" x ")[1]) * u.s

    def _read_obsid(self, obsid):
        """
        Reads different fields from OBS ID number.
        """
        options = {}
        data = {"obsid": obsid}
        field_keys = {
            "Large linelist": "linelist",
            "Default compression": "compression",
            "Lossy compression": "compression",
            "Non-simultaneous readout": "readout",
            "SJI cadence default": "sji_cadence",
            "SJI cadence 10s": "sji_cadence",
            "FUV binned same as NUV": "fuv_binning",
            "Spatial x 1, Spectral x 1": "binning",
            "Exposure 1s": "exptime",
            "C II   Si IV   Mg II h/k   Mg II w": "sjis",
        }
        versions = [34, 36, 38, 40]
        if len(str(obsid)) != 10:
            msg = "Invalid OBS ID: must have 10 digits."
            raise ValueError(msg)
        # iris_obsid2struct.pro: the table version is obsid // 2e8 * 2, so v40 Full readout OBS IDs start with 41
        version = obsid // 200_000_000 * 2
        if version >= 42:
            msg = (
                f"OBS ID {obsid} is a calibration or special sequence, which has no OBS table. "
                "The OBS_DESC and OBSTITLE keywords in the file header describe it."
            )
            raise ValueError(msg)
        if version not in versions:
            msg = f"Invalid OBS ID: OBS table version {version} is not one of {versions}"
            raise ValueError(msg)
        obsid -= version * 100_000_000
        reader = ascii_io.get_reader(reader_cls=ascii_io.Csv)
        with importlib.resources.as_file(importlib.resources.files("irispy") / f"data/v{version}-table10.csv") as path:
            table1 = reader.read(path)
        with importlib.resources.as_file(
            importlib.resources.files("irispy") / f"data/v{version}-table2000.csv"
        ) as path:
            table2 = reader.read(path)
        id_raster = int(str(obsid)[-2:])
        try:
            meta = table1[table1["OBS-ID"] == id_raster][0]
        except IndexError:
            msg = f"Invalid OBS ID: last two numbers must be between {table1['OBS-ID'].min()} and {table1['OBS-ID'].max()}"
            raise ValueError(
                msg,
            ) from None

        # The v40 sparse/dense/sparse rasters give one step per sub-raster, "1  0.33  1".
        data["raster_step"] = np.squeeze(np.array(str(meta["Raster step"]).split(), dtype=float)) * u.arcsec
        # The v40 table writes "126 x120" where the others write "126x120".
        data["raster_fov"] = str(meta["Raster FOV"]).replace(" x", "x")
        data["spec_cadence"] = u.Quantity(meta["Spectral cadence"])
        data["sji_fov"] = str(meta["SJI FOV"]).replace(" x", "x")
        data["raster_desc"] = str(meta["Description"])
        data["raster_fulldesc"] = f"{data['raster_desc']} {data['raster_fov']} {meta['Spectral cadence']}"
        field_ranges = np.r_[np.flatnonzero(table2["OBS ID"] == 0), len(table2)]
        # field indices, start from largest and subtract
        for start, end in zip(field_ranges[-2::-1], field_ranges[:0:-1], strict=True):
            table = table2[start:end]
            for i in np.arange(start, end)[::-1]:
                index = i
                tmp = table2["OBS ID"][i]
                if (obsid - tmp) >= 0:
                    obsid -= tmp
                    break
            desc = table["Size + description"].tolist()
            # Save values for attributes but also table options as function of OBS ID
            if desc[0] in field_keys:
                attr_name = field_keys[desc[0]]
                if attr_name == "exptime":
                    opt = (self._exptime_to_quant(a) for a in desc)
                    opt = dict(zip(opt, table["OBS ID"].tolist(), strict=True))
                    attr_value = self._exptime_to_quant(desc[index - start])
                else:
                    opt = dict(zip(desc, table["OBS ID"].tolist(), strict=True))
                    attr_value = desc[index - start]
                data[attr_name] = attr_value
                options[attr_name] = opt
        if obsid != id_raster:
            msg = f"Invalid OBS ID: no combination of v{version} OBS table options gives {data['obsid']}"
            raise ValueError(msg)
        return data, options
