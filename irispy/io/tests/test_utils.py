from irispy.io.utils import fits_info, read_files, tar_extract_all
from irispy.sji import SJICube
from irispy.spectrograph import Collection


def test_tar_extract_all(fake_tar):
    # fake_tar is a Path
    extracted_files = tar_extract_all(fake_tar)
    assert len(extracted_files) == 4
    for file in extracted_files:
        assert file.exists()
        assert file.is_file()
    assert all(file.suffix == ".fits" for file in extracted_files)
    assert all(file.parent == fake_tar.parent / fake_tar.stem for file in extracted_files)


def test_fits_info(capsys, raster_file, sji_1330_file, sji_1400_file, sji_2796_file, sji_2832_file):
    fits_info(raster_file)
    captured = capsys.readouterr()
    assert raster_file in captured.out

    fits_info(sji_1330_file)
    captured = capsys.readouterr()
    assert sji_1330_file in captured.out

    fits_info(sji_1400_file)
    captured = capsys.readouterr()
    assert sji_1400_file in captured.out

    fits_info(sji_2796_file)
    captured = capsys.readouterr()
    assert sji_2796_file in captured.out

    fits_info(sji_2832_file)
    captured = capsys.readouterr()
    assert sji_2832_file in captured.out


def test_read_files_raster(raster_file):
    output = read_files(raster_file)
    assert output["raster"] is not None
    assert len(output["raster"]) == 1
    assert not output["sji"]

    output = read_files([raster_file])
    assert output["raster"] is not None
    assert len(output["raster"]) == 1
    assert not output["sji"]

    output = read_files([raster_file, raster_file])
    assert output["raster"] is not None
    assert len(output["raster"]) == 2
    assert not output["sji"]


def test_read_files_sji(sji_1330_file, sji_1400_file, sji_2796_file, sji_2832_file):
    output = read_files(sji_1330_file)
    assert output["sji"] is not None
    assert len(output["sji"]) == 1
    assert not output["raster"]

    output = read_files(sji_1400_file)
    assert output["sji"] is not None
    assert len(output["sji"]) == 1
    assert not output["raster"]

    output = read_files(sji_2796_file)
    assert output["sji"] is not None
    assert len(output["sji"]) == 1
    assert not output["raster"]

    output = read_files(sji_2832_file)
    assert output["sji"] is not None
    assert len(output["sji"]) == 1
    assert not output["raster"]

    output = read_files([sji_2832_file, sji_2796_file, sji_1400_file, sji_1330_file])
    assert output["sji"] is not None
    assert len(output["sji"]) == 4
    assert not output["raster"]


def test_read_files_with_mix(raster_file, sji_1330_file):
    output = read_files([raster_file, sji_1330_file])
    assert output["sji"] is not None
    assert len(output["sji"]) == 1
    assert isinstance(output["sji"][0], SJICube)

    assert output["raster"] is not None
    assert len(output["raster"]) == 1
    assert isinstance(output["raster"][0], Collection)
