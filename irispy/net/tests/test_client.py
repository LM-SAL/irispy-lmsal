import json
from pathlib import Path
from unittest import mock

import numpy as np
import parfive
import pytest
import requests

import astropy.units as u

from sunpy.net import Fido
from sunpy.net import attrs as a
from sunpy.util.datatype_factory_base import NoMatchError

from irispy.data.test import get_test_filepath
from irispy.io import read_files
from irispy.net import LMSALClient, client

TIME = a.Time("2020-03-03 00:00", "2020-03-03 05:00")


@pytest.fixture
def hcr(monkeypatch):
    """
    Answer the HCR searches with a recorded response, returning the mock of
    `requests.get`.
    """

    def record(name):
        response = mock.Mock()
        response.json.return_value = json.loads(Path(get_test_filepath(f"hcr/{name}")).read_text(encoding="latin-1"))
        get = mock.Mock(return_value=response)
        monkeypatch.setattr(requests, "get", get)
        return get

    return record


def test_search(hcr):
    get = hcr("hcr_20200303.json")
    table = LMSALClient().search(TIME)
    assert get.call_args.kwargs["params"] == {
        "cmd": "search-events3",
        "outputformat": "json",
        "instrument": "IRIS",
        "hasData": "true",
        "limit": 5000,
        "startTime": "2020-03-03T00:00:00.000",
        "stopTime": "2020-03-03T05:00:00.000",
    }
    # The rasters and slit-jaw files of the three observations, in time order; not the mosaic.
    assert list(table["Filename"]) == [
        "iris_l2_20200303_010044_3640604076_SJI_1330_t000.fits.gz",
        "iris_l2_20200303_010044_3640604076_SJI_1400_t000.fits.gz",
        "iris_l2_20200303_010044_3640604076_SJI_2796_t000.fits.gz",
        "iris_l2_20200303_010044_3640604076_SJI_2832_t000.fits.gz",
        "iris_l2_20200303_010044_3640604076_raster.tar.gz",
        "iris_l2_20200303_041957_3882010194_SJI_1330_t000.fits.gz",
        "iris_l2_20200303_041957_3882010194_SJI_1400_t000.fits.gz",
        "iris_l2_20200303_041957_3882010194_SJI_2796_t000.fits.gz",
        "iris_l2_20200303_041957_3882010194_SJI_2832_t000.fits.gz",
        "iris_l2_20200303_041957_3882010194_raster.tar.gz",
        "iris_l2_20200303_044800_3605202860_SJI_1400_t000.fits.gz",
        "iris_l2_20200303_044800_3605202860_SJI_2796_t000.fits.gz",
        "iris_l2_20200303_044800_3605202860_SJI_2832_t000.fits.gz",
        "iris_l2_20200303_044800_3605202860_raster.tar.gz",
    ]
    raster, sji = table[-1], table[-4]
    assert raster["OBSID"] == 3605202860
    assert u.allclose([raster["FOV X"], raster["FOV Y"]], [63.4, 175.1] * u.arcsec)
    assert np.ma.is_masked(raster["Wavelength"])
    assert sji["Wavelength"] == 1400 * u.AA
    # A slit-jaw file has its own times; the observation started earlier
    assert sji["Start Time"].isot == "2020-03-03T04:48:06.000"
    assert sji["Observation Start"].isot == raster["Start Time"].isot == "2020-03-03T04:48:00.000"
    header = str(table).splitlines()[0]
    assert "Filename" in header
    assert "Cadence" not in header
    # No observation has a NOAA active region number
    assert table["NOAA AR"].mask.all()


def test_search_attrs(hcr):
    get = hcr("hcr_20200303.json")
    Fido.search(
        TIME,
        a.Instrument.iris,
        a.iris.ObsID("3605202860"),
        a.iris.Target("qs"),
        a.iris.XCen(-1 * u.arcmin, 0 * u.arcmin),
        a.iris.YCen(0 * u.arcsec, 20 * u.arcsec),
    )
    params = get.call_args.kwargs["params"]
    assert params["obsid"] == 3605202860
    assert params["target"] == "QS"
    assert (params["minxcen"], params["maxxcen"]) == (-60, 0)
    assert (params["minycen"], params["maxycen"]) == (0, 20)


@pytest.mark.parametrize(
    ("text", "observations"),
    [
        # In the title
        ("Throughput", {"20200303_041957_3882010194"}),
        # In the goal, the OBS description
        ("DENSE", {"20200303_010044_3640604076"}),
    ],
)
def test_search_description(hcr, text, observations):
    hcr("hcr_20200303.json")
    table = Fido.search(TIME, a.Provider("LMSAL"), a.iris.Description(text))["lmsal"]
    assert {name[8:34] for name in table["Filename"]} == observations


def test_search_products(hcr):
    hcr("hcr_20160108.json")
    time = a.Time("2016-01-08 18:00", "2016-01-08 20:00")
    sot = Fido.search(time, a.Provider("LMSAL"), a.iris.Product("SOT"))["lmsal"]
    assert list(sot["Filename"]) == [
        "iris_l2_20160108_173544_3680100932_SOTFG.tar.gz",
        "iris_l2_20160108_191211_3680100932_SOTFG.tar.gz",
        "iris_l2_20160108_191211_3680100932_SOTSP.tar.gz",
    ]
    sdo = Fido.search(time, a.Provider("LMSAL"), a.iris.Product.sdo)["lmsal"]
    assert list(sdo["Filename"]) == [
        "iris_l2_20160108_173544_3680100932_SDO.tar.gz",
        "iris_l2_20160108_191211_3680100932_SDO.tar.gz",
    ]
    # The tarballs have the observation's times and no field of view or wavelength.
    assert [time.isot for time in sdo["Start Time"]] == ["2016-01-08T17:35:44.000", "2016-01-08T19:12:11.000"]
    for name in ["FOV X", "FOV Y", "Wavelength", "Images", "Cadence"]:
        assert sdo[name].mask.all()
    assert sot.colnames == sdo.colnames


def test_search_deconvolved(hcr):
    hcr("hcr_20200303.json")
    table = Fido.search(TIME, a.Provider("LMSAL"), a.iris.Product("deconvolved"))["lmsal"]
    assert list(table["Filename"]) == [
        "iris_l2_20200303_041957_3882010194_SJI_2796_t000_deconvolved.fits.gz",
        "iris_l2_20200303_041957_3882010194_SJI_2832_t000_deconvolved.fits.gz",
    ]
    assert list(table["Product"]) == ["deconvolved", "deconvolved"]


def test_search_limit(hcr, monkeypatch):
    hcr("hcr_20200303.json")
    monkeypatch.setattr(client, "HCR_LIMIT", 4)
    with pytest.warns(UserWarning, match="returns at most 4 observations"):
        LMSALClient().search(TIME)


@pytest.mark.parametrize(
    "query",
    [
        (a.Instrument.iris,),
        (a.Provider("LMSAL"),),
        (a.Provider("lmsal"),),
    ],
)
def test_search_instrument_provider(hcr, monkeypatch, query):
    hcr("hcr_20200303.json")
    # Only LMSALClient, so that no other client searches.
    monkeypatch.setattr(Fido, "registry", {LMSALClient: LMSALClient._can_handle_query})
    result = Fido.search(TIME, *query)
    assert result.keys() == ["lmsal"]
    assert len(result["lmsal"]) == 14


@pytest.mark.parametrize("query", [(TIME,), (TIME, a.Instrument.aia), (TIME, a.Provider("VSO"))])
def test_search_not_iris(hcr, monkeypatch, query):
    get = hcr("hcr_20200303.json")
    # Only LMSALClient, so that no other client searches.
    monkeypatch.setattr(Fido, "registry", {LMSALClient: LMSALClient._can_handle_query})
    with pytest.raises(NoMatchError):
        Fido.search(*query)
    get.assert_not_called()


def test_search_iris_attr_only(hcr):
    get = hcr("hcr_20200303.json")
    # No client answers: the a.iris attrs only refine a search for IRIS, and the VSO does not know them.
    with pytest.raises(NoMatchError):
        Fido.search(TIME, a.iris.Target("AR"))
    get.assert_not_called()


def test_fetch(hcr, tmp_path):
    hcr("hcr_20200303.json")
    result = Fido.search(TIME, a.Provider("LMSAL"), a.iris.Product("raster"))
    downloader = mock.create_autospec(parfive.Downloader, instance=True)
    Fido.fetch(result, path=tmp_path, downloader=downloader)
    assert [(call.args[0], call.kwargs["filename"]) for call in downloader.enqueue_file.call_args_list] == [
        (
            f"https://www.lmsal.com/solarsoft/irisa/data/level2_compressed/2020/03/03/{obs}/iris_l2_{obs}_raster.tar.gz",
            str(tmp_path / obs / f"iris_l2_{obs}_raster.tar.gz"),
        )
        for obs in ["20200303_010044_3640604076", "20200303_041957_3882010194", "20200303_044800_3605202860"]
    ]


@pytest.mark.remote_data
def test_search_fetch_remote(tmp_path):
    result = Fido.search(a.Time("2021-10-01 06:09", "2021-10-01 06:12"), a.Instrument.iris, a.iris.Product.sji)
    table = result["lmsal"]
    assert set(table["OBSID"]) == {3683602040}
    files = Fido.fetch(
        table[table["Filename"] == "iris_l2_20211001_060925_3683602040_SJI_1400_t000.fits.gz"], path=tmp_path
    )
    assert files == [
        str(tmp_path / "20211001_060925_3683602040" / "iris_l2_20211001_060925_3683602040_SJI_1400_t000.fits.gz")
    ]
    assert read_files(tmp_path / "20211001_060925_3683602040").meta.observing_mode_id == 3683602040
