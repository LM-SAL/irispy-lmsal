import re
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import animation, patheffects

from astropy.visualization import AsinhStretch, AsymmetricPercentileInterval, ImageNormalize

from sunpy.visualization.colormaps.color_tables import iris_sji_color_table

from .utils import _boolean_mask, _saturated_samples

__all__ = ["generate_wobble_movie"]


def generate_wobble_movie(
    cube,
    *,
    outdir: str | Path = "./",
    trim: bool = False,
    timestamp: bool = True,
    wobble_cadence: int = 180,
    ffmpeg_path: str | Path | None = None,
    **kwargs,
) -> Path:
    """
    Creates a wobble movie from a slit-jaw cube.

    ..note:

        This requires FFMPEG to be installed and discoverable.
        If FFMPEG is not found, you can pass it as an argument called ``ffmpeg_path``.

    Parameters
    ----------
    cube : `~irispy.sji.SJICube`
        The slit-jaw images to create a wobble movie from.
    outdir : Union[str, pathlib.Path], optional
        Location to save the movie.
        Defaults to the current working directory.
    trim : `bool`, optional
        Movie is trimmed to the rows and columns that have data in all frames, by default False.
    timestamp : `bool`, optional
        If `True`, will add the time of each frame to the wobble movie.
        Optional, defaults to `True`.
    wobble_cadence : `int`, optional
        Sets the cadence of the wobble movie in seconds.
        The movie shows every n-th frame of the cube, n being ``wobble_cadence`` divided by the median
        time between frames, rounded down and at least 1.
        Optional, defaults to 180 seconds.
    ffmpeg_path : Union[str, pathlib.Path], optional
        Path to FFMPEG executable, by default `None`.
        In theory you will not need to do this but matplotlib might not be able to find the ffmpeg exe.
    **kwargs : `dict`, optional
        Keyword arguments to passed to `matplotlib.animation.FuncAnimation`.

    Returns
    -------
    `pathlib.Path`
        The movie created.

    Raises
    ------
    ValueError
        If ``trim`` is `True` and no pixel has data in every frame, as in observations with a moving
        field of view.

    Notes
    -----
    - 2832 is considered the best wavelength to use for wobble movies.
    """
    if ffmpeg_path:
        mpl.rcParams["animation.ffmpeg_path"] = ffmpeg_path
    # Trim down to the rows and columns that contain data in all frames
    if trim:
        # The reader masks saturated samples, which still hold data
        has_data = ~(_boolean_mask(cube) & ~_saturated_samples(cube)).any(axis=0)
        if not has_data.any():
            msg = "No pixel has data in every frame, so there is nothing to trim to; use trim=False."
            raise ValueError(msg)
        rows, columns = has_data.sum(axis=1), has_data.sum(axis=0)
        (rows,) = np.where(rows > rows.max() * 0.8)
        (columns,) = np.where(columns > columns.max() * 0.8)
        cube = cube[:, rows[0] : rows[-1] + 1, columns[0] : columns[-1] + 1]
    times = cube.time
    # Downsample in time to accentuate the wobble
    step = 1
    if len(times) > 1:
        cadence = np.median((times[1:] - times[:-1]).to_value("s"))
        step = max(1, int(wobble_cadence // cadence))
    frames = range(0, len(times), step)

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection=cube[0].wcs)
    colormap = iris_sji_color_table(str(int(cube.meta["TWAVE1"])))
    # 0.15% at each end, the default cutoff of the image_clipping irispy used here before (V. Hansteen's
    # port of an IDL routine by P. Suetterlin)
    vmin, vmax = AsymmetricPercentileInterval(0.15, 99.85).get_limits(cube.data[~_boolean_mask(cube)])
    image = ax.imshow(
        cube.data[0],
        origin="lower",
        cmap=colormap,
        norm=ImageNormalize(vmin=vmin, vmax=vmax, stretch=AsinhStretch()),
    )
    ax.set_xlabel("Solar X")
    ax.set_ylabel("Solar Y")
    if timestamp:
        title = ax.text(
            0.5,
            0.95,
            str(times[0]),
            color="w",
            transform=ax.transAxes,
            ha="center",
            path_effects=[patheffects.withStroke(linewidth=3, foreground="black")],
        )
    else:
        title = ax.text(0.5, 0.95, "")

    def update(i):
        image.set_array(cube.data[i])
        if timestamp:
            title.set_text(str(times[i]))
        return image, title

    anim = animation.FuncAnimation(fig, func=update, frames=frames, blit=True, repeat=False, **kwargs)
    clean_filename = re.sub(r"[^\w\-_\. ]", "_", f"{cube.meta['TDESC1']}_{times[0].isot.split('.')[0]}_wobble.mp4")
    filename = Path(outdir) / clean_filename
    anim.save(filename, writer=animation.FFMpegWriter(fps=12))
    return filename
