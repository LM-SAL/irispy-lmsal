from typing import Union, Optional
from pathlib import Path

import matplotlib.animation as animation
import matplotlib.patheffects as PathEffects
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.time import TimeDelta
from astropy.visualization import AsinhStretch, ImageNormalize
from astropy.wcs import WCS
from sunpy.time import parse_time
from sunpy.visualization.colormaps.color_tables import iris_sji_color_table

from irispy.utils import image_clipping

__all__ = ["wobble_movie"]


def wobble_movie(
    filelist: list,
    outdir: Union[str, Path] = "./",
    trim: bool = False,
    timestamp: bool = False,
    wobble_cadence: int = 180,
    ffmpeg_path: Optional[Union[str, Path]] = None,
    **kwargs,
) -> None:
    """
    Creates a wobble movie from a list of files.

    ..note:

        This requires FFMPEG to be installed and discoverable.
        If FFMPEG is not found, you can pass it as an argument called ``ffmpeg_path``.

    Parameters
    ----------
    filelist : `list`
        Files to create a wobble movie from.
    outdir : Union[str,Path], optional
        Location to save the movie(s).
        Defaults to the current working directory.
    trim : `bool`, optional
        Movie is trimmed to include only area that has data in all frames, by default False
    timestamp : `bool`, optional
        If `True`, will add a timestamp to the wobble movie.
        Optional, defaults to `False`.
    wobble_cadence : `int`, optional
        Sets the cadence of the wobble movie in seconds.
        Optional, defaults to 180 seconds.
    ffmpeg_path : Union[str,Path], optional
        Path to FFMPEG executable, by default `None`.
        In theory you will not need to do this but matplotlib might not be able to find the ffmpeg exe.
    **kwargs : `dict`, optional
        Keyword arguments to passed to `FuncAnimation`.

    Returns
    -------
    `list`
        A list of the movies created.

    Notes
    -----
    This is designed to be used on IRIS Level 2 SJI data.

    2832 is considered the best wavelength to use for wobble movies.

    Timestamps take the main header cadence and add that to the "DATEOBS".
    They do not use the information in the AUX array.
    """
    if ffmpeg_path:
        import matplotlib as mpl

        mpl.rcParams["animation.ffmpeg_path"] = ffmpeg_path

    filenames = []
    for file in filelist:
        data, header = fits.getdata(file, header=True)
        wcs = WCS(header)
        numframes = header["NAXIS3"]
        date = header["DATE_OBS"].split(".")[0]
        # Calculate index to downsample in time to accentuate the wobble
        cadence = header["CDELT3"]
        cadence_sample = np.floor(wobble_cadence / cadence) if np.floor(wobble_cadence / cadence) > 1 else 1
        if timestamp:
            timestamps = [
                parse_time(header["STARTOBS"]) + TimeDelta(cadence, format="sec") * i for i in range(numframes)
            ]
        else:
            timestamps = [parse_time(header["STARTOBS"])]
        # Trim down to only that part of the movie that contains data in all frames
        if trim:
            # TODO: improve this, it trims a bit but not fully
            dmin = np.min(data, axis=0)
            dmask = dmin > -200
            dmx = np.sum(dmask, axis=1)
            dmy = np.sum(dmask, axis=0)
            (subx,) = np.where(dmx > (np.max(dmx) * 0.8))
            (suby,) = np.where(dmy > (np.max(dmy) * 0.8))
            data = data[:, suby[0] : suby[-1], subx[0] : subx[-1]]

        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1, projection=wcs.dropaxis(-1))
        colormap = iris_sji_color_table(str(int(header["TWAVE1"])))
        vmin, vmax = image_clipping(data)
        image = ax.imshow(
            data[0],
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
                str(timestamps[0]),
                color="w",
                transform=ax.transAxes,
                ha="center",
                path_effects=[PathEffects.withStroke(linewidth=3, foreground="black")],
            )
        else:
            title = ax.text(0.5, 0.95, "")

        def update(i):
            image.set_array(data[i])
            if timestamp:
                title.set_text(str(timestamps[i]))
            return image, title

        anim = animation.FuncAnimation(
            fig, func=update, frames=range(0, numframes, int(cadence_sample)), blit=True, repeat=False, **kwargs
        )
        filename = Path(outdir) / Path(f"{header['TDESC1']}_{date}_wobble.mp4")
        writervideo = animation.FFMpegWriter(fps=12)
        anim.save(filename, writer=writervideo)
        filenames.append(filename)
    return filenames
