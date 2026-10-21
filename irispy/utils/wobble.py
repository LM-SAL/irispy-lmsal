import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import animation, patheffects

import astropy.units as u
from astropy.visualization import AsinhStretch, AsymmetricPercentileInterval, ImageNormalize

from sunpy.visualization.colormaps.color_tables import iris_sji_color_table

from .utils import _boolean_mask, _saturated_samples

__all__ = ["generate_wobble_movie"]


@u.quantity_input
def generate_wobble_movie(
    cube,
    *,
    outdir: str | Path = "./",
    trim: bool = False,
    timestamp: bool = True,
    wobble_cadence: u.s = 180 * u.s,
) -> Path:
    """
    Creates a wobble movie from a slit-jaw cube.

    .. note::

        This requires FFmpeg. If matplotlib does not find it, set
        ``matplotlib.rcParams["animation.ffmpeg_path"]`` to the path of the FFmpeg executable.

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
    wobble_cadence : `astropy.units.Quantity`, optional
        Sets the cadence of the wobble movie.
        The movie shows every n-th frame of the cube, n being ``wobble_cadence`` divided by the median
        time between frames, rounded to the nearest integer (halves up) and at least 1.
        Optional, defaults to 180 seconds.

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
    `iris_check_wobble.pro <https://sohoftp.nascom.nasa.gov/solarsoft/iris/idl/ucar/iris_check_wobble.pro>`__
    recommends the 2832 passband for wobble movies.
    """
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
        # Halves round up, as IDL's ROUND does in iris_check_wobble.pro (Python's round gives round(2.5) == 2)
        step = max(1, int(wobble_cadence.to_value(u.s) / cadence + 0.5))
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

    anim = animation.FuncAnimation(fig, func=update, frames=frames, blit=True, repeat=False)
    clean_filename = re.sub(r"[^\w\-_\. ]", "_", f"{cube.meta['TDESC1']}_{times[0].isot.split('.')[0]}_wobble.mp4")
    filename = Path(outdir) / clean_filename
    try:
        anim.save(filename, writer=animation.FFMpegWriter(fps=12))
    finally:
        plt.close(fig)
    return filename
