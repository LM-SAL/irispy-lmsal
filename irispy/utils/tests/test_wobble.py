import shutil
import subprocess

import matplotlib.pyplot as plt
import numpy as np
import pytest
from matplotlib import animation

import astropy.units as u

from irispy.utils.wobble import generate_wobble_movie


@pytest.fixture
def figures(monkeypatch):
    # The figures generate_wobble_movie draws, which it closes once the movie is saved
    figures = []
    figure = plt.figure
    monkeypatch.setattr(plt, "figure", lambda: figures.append(figure()) or figures[-1])
    return figures


@pytest.mark.parametrize(("trim", "spatial"), [(False, np.s_[:, :]), (True, np.s_[1:38, 1:36])])
def test_generate_wobble_movie(sns_sjicube_2832, tmp_path, figures, trim, spatial):
    # The frames are about 1900 s apart, so a 1.1 hour (3960 s) cadence shows every second frame
    movie = generate_wobble_movie(sns_sjicube_2832, outdir=tmp_path, trim=trim, wobble_cadence=1.1 * u.hour)

    assert movie.parent == tmp_path
    ffmpeg = shutil.which(animation.FFMpegWriter.bin_path())
    decoded = subprocess.run(  # noqa: S603 - decode a movie created in this test's temporary directory.
        [
            ffmpeg,
            "-v",
            "error",
            "-i",
            str(movie),
            "-vf",
            "scale=1:1",
            "-pix_fmt",
            "gray",
            "-f",
            "rawvideo",
            "-",
        ],
        capture_output=True,
        check=True,
    )
    # Each decoded frame is one grayscale pixel, so one byte represents one frame.
    assert len(decoded.stdout) == 5
    (figure,) = figures
    ax = figure.axes[0]
    assert ax.texts[0].get_text() == sns_sjicube_2832.time[8].isot
    np.testing.assert_array_equal(ax.images[0].get_array(), sns_sjicube_2832[(8, *spatial)].data)
    # Solar X and Y follow the trim
    np.testing.assert_allclose(
        ax.wcs.pixel_to_world_values(0, 0),
        sns_sjicube_2832[(0, *spatial)].wcs.low_level_wcs.pixel_to_world_values(0, 0),
    )


def test_generate_wobble_movie_one_frame(sns_sjicube_2832, tmp_path, figures):
    assert generate_wobble_movie(sns_sjicube_2832[:1], outdir=tmp_path).exists()
    # The figure is closed, so making more than 20 movies does not trip matplotlib's open-figure warning
    assert not plt.get_fignums()
    np.testing.assert_array_equal(figures[0].axes[0].images[0].get_array(), sns_sjicube_2832[0].data)


def test_generate_wobble_movie_trim_without_common_data(sns_sjicube_2832, tmp_path):
    # A frame without data leaves no pixel with data in every frame
    sns_sjicube_2832.mask[0] = True
    with pytest.raises(ValueError, match="No pixel has data in every frame"):
        generate_wobble_movie(sns_sjicube_2832, outdir=tmp_path, trim=True)
