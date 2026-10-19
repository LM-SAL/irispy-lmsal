import sunpy

from irispy.data import sample


def test_file_dict_finds_samples_in_sunpy_sample_dir(tmp_path, monkeypatch):
    monkeypatch.delenv("SUNPY_SAMPLEDIR", raising=False)
    monkeypatch.setitem(sunpy.config["downloads"], "sample_dir", str(tmp_path))
    sji_1330 = tmp_path / "iris_l2_20211001_060925_3683602040_SJI_1330_t000.fits.gz"
    sji_1330.touch()
    assert sample.file_dict == {
        "RASTER_FITS": None,
        "RASTER_TAR": None,
        "SJI_1330": sji_1330,
        "SJI_1400": None,
        "SJI_2796": None,
        "SJI_2832": None,
    }
    assert sample.file_list == [sji_1330]
