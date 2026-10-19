"""
Convert the IRIS pipeline's bad-pixel tables into
``irispy/data/iris_sji_badpix_<date>.fits``.

The pipeline keeps the dust and dead pixels of each slit-jaw channel as a list
of detector addresses in ``<stamp>_badpix.geny`` (an IDL save file) indexed by
the flat-field table ``<stamp>_flat.genx``; both live in the SolarSoft database
``$SSWDB/iris/data`` with a public copy at ``SRC_URL``. ``IRIS_PREP_GET_BADPIX``
picks the record nearest the observation time, but every record of a channel
has held the same pixel list since 2014, so the table is written once per
channel and this script stops if a new file breaks that assumption.

Run by hand when a new pair of files appears on the server:

    python tools/badpix.py 20260822_041540
"""

import argparse
import datetime
import tempfile
import urllib.request
from pathlib import Path

import numpy as np
from scipy.io import readsav

from astropy.io import fits

from sunpy.io.special.genx import read_genx

SRC_URL = "https://sohoftp.nascom.nasa.gov/sdb/iris/data/"
DETECTOR_COLUMNS, DETECTOR_ROWS = 2072, 1096
# iris_mk_pointdb.pro (version 11): the slit center on each channel's detector
# (1-based FITS pixel), the plate scale (arcsec per pixel) and the angle of the
# detector y axis east of the slit (degrees).
POINTING = {
    "SJI_1330": (537.30, 524.03, 0.16560, -0.200),
    "SJI_1400": (529.32, 510.46, 0.16560, -0.224),
    "SJI_1600W": (534.0, 516.7, 0.16560, -0.21),
    "SJI_2796": (504.69, 503.402, 0.16790, 0.274),
    "SJI_2832": (506.47, 502.22, 0.16790, 0.286),
    "SJI_5000W": (505.2, 504.7, 0.16790, 0.28),
}
# cdlt_1p5 in iris_mk_pointdb.pro: the level 2 plate scale of an unsummed pixel.
LEVEL2_CDELT = 0.16635


def fetch(stamp, directory):
    paths = []
    for suffix in ("badpix.geny", "flat.genx"):
        path = Path(directory) / f"{stamp}_{suffix}"
        if not path.exists():
            print(f"downloading {SRC_URL}{path.name}")
            urllib.request.urlretrieve(SRC_URL + path.name, path)  # NOQA: S310
        paths.append(path)
    return paths


def pixel_sets(badpix_path, flat_path):
    records = readsav(str(badpix_path))["p0"][0]
    sets = {}
    for row in read_genx(str(flat_path))["SAVEGEN0"]:
        channel = str(row["IMG_PATH"])
        if channel not in POINTING:
            continue
        addresses = np.unique(np.asarray(records[f"F{int(row['RECNUM'])}"], dtype=np.int64))
        if not np.array_equal(sets.setdefault(channel, addresses), addresses):
            msg = f"{channel}: record {int(row['RECNUM'])} differs from the earlier records; the table is no longer one list per channel"
            raise SystemExit(msg)
    return sets


def _pixel_hdu(channel, addresses):
    x, y = addresses % DETECTOR_COLUMNS, addresses // DETECTOR_COLUMNS
    order = np.lexsort((y, x))
    hdu = fits.BinTableHDU.from_columns(
        [
            fits.Column("X", "I", array=x[order].astype(np.int16)),
            fits.Column("Y", "I", array=y[order].astype(np.int16)),
        ],
        name=channel,
    )
    hdu.header.comments["TTYPE1"] = "detector column, 0-based, before summing"
    hdu.header.comments["TTYPE2"] = "detector row, 0-based, before summing"
    return hdu


def write(sets, output, badpix_name, flat_name):
    channels = sorted(sets)
    pointing = fits.BinTableHDU.from_columns(
        [
            fits.Column("CHANNEL", "9A", array=channels),
            fits.Column("CPX1", "D", array=[POINTING[c][0] for c in channels]),
            fits.Column("CPX2", "D", array=[POINTING[c][1] for c in channels]),
            fits.Column("CDELT", "D", array=[POINTING[c][2] for c in channels]),
            fits.Column("ROLL", "D", array=[POINTING[c][3] for c in channels]),
        ],
        name="POINTING",
    )
    pointing.header.comments["TTYPE2"] = "slit center column, 1-based FITS pixel"
    pointing.header.comments["TTYPE3"] = "slit center row, 1-based FITS pixel"
    pointing.header.comments["TTYPE4"] = "plate scale, arcsec per detector pixel"
    pointing.header.comments["TTYPE5"] = "detector y axis east of the slit, degrees"
    primary = fits.PrimaryHDU()
    primary.header["SRCBADPX"] = (badpix_name, "pipeline bad-pixel file in $SSWDB/iris/data")
    primary.header["SRCFLAT"] = (flat_name, "flat-field index keying the bad-pixel file")
    primary.header["SRCURL"] = (SRC_URL, "public copy of $SSWDB/iris/data")
    primary.header["DETNX"] = (DETECTOR_COLUMNS, "detector columns")
    primary.header["DETNY"] = (DETECTOR_ROWS, "detector rows")
    primary.header["L2CDELT"] = (LEVEL2_CDELT, "level 2 plate scale, arcsec per unsummed pixel")
    primary.header["PDBVER"] = (11, "iris_mk_pointdb.pro version of POINTING")
    primary.header["DATE"] = (datetime.datetime.now(tz=datetime.UTC).date().isoformat(), "conversion date")
    pixels = [_pixel_hdu(c, sets[c]) for c in channels]
    fits.HDUList([primary, pointing, *pixels]).writeto(output, overwrite=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("stamp", help="file stamp on the server, for example 20260822_041540")
    parser.add_argument(
        "--source", help="directory holding the two files; downloaded to a temporary directory otherwise"
    )
    args = parser.parse_args()
    output = Path(__file__).resolve().parents[1] / "irispy" / "data" / f"iris_sji_badpix_{args.stamp[:8]}.fits"
    with tempfile.TemporaryDirectory() as tmp:
        badpix_path, flat_path = fetch(args.stamp, args.source or tmp)
        sets = pixel_sets(badpix_path, flat_path)
    write(sets, output, badpix_path.name, flat_path.name)
    for channel in sorted(sets):
        print(f"{channel}: {len(sets[channel])} pixels")
    print(f"wrote {output} ({output.stat().st_size / 1e3:.0f} kB)")


if __name__ == "__main__":
    main()
