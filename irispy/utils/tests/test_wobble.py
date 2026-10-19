import shutil
import subprocess

import matplotlib.pyplot as plt
import numpy as np
import pytest
from matplotlib import animation

from irispy.utils.wobble import generate_wobble_movie


@pytest.mark.parametrize(("trim", "spatial"), [(False, np.s_[:, :]), (True, np.s_[1:38, 1:36])])
def test_generate_wobble_movie(sns_sjicube_2832, tmp_path, trim, spatial):
    # The frames are about 1900 s apart, so a 4000 s cadence shows every second frame
    movie = generate_wobble_movie(sns_sjicube_2832, outdir=tmp_path, trim=trim, wobble_cadence=4000)

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
    ax = plt.gcf().axes[0]
    assert ax.texts[0].get_text() == sns_sjicube_2832.time[8].isot
    np.testing.assert_array_equal(ax.images[0].get_array(), sns_sjicube_2832[(8, *spatial)].data)
    # Solar X and Y follow the trim
    np.testing.assert_allclose(
        ax.wcs.pixel_to_world_values(0, 0),
        sns_sjicube_2832[(0, *spatial)].wcs.low_level_wcs.pixel_to_world_values(0, 0),
    )
    plt.close("all")


def test_generate_wobble_movie_trim_without_common_data(sns_sjicube_2832, tmp_path):
    # A frame without data leaves no pixel with data in every frame
    sns_sjicube_2832.mask[0] = True
    with pytest.raises(ValueError, match="No pixel has data in every frame"):
        generate_wobble_movie(sns_sjicube_2832, outdir=tmp_path, trim=True)
