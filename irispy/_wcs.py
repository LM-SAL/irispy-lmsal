import numpy as np

from astropy.utils.shapes import simplify_basic_index
from astropy.wcs import WCS
from astropy.wcs.utils import wcs_to_celestial_frame
from astropy.wcs.wcsapi.wrappers.sliced_wcs import sanitize_slices

from sunpy.coordinates.frames import Helioprojective


def _celestial_frame_from_cube(cube):
    """
    The `~sunpy.coordinates.frames.Helioprojective` frame of this observation, as seen
    by the IRIS observer.

    This does not depend on slicing: it works on combined multi-file cubes and on
    cubes sliced along the scan, step, slit, time, or wavelength axes.
    """
    observer = getattr(cube.meta, "observer", None)
    if observer is not None:
        return Helioprojective(observer=observer, obstime=observer.obstime)

    headers = cube.meta.get("frame_wcs_headers")
    if headers is not None:
        # The first frame, whatever the slice: fits_wcs follows an integer index and a rebinned cube has none.
        fits_wcs = WCS(headers[0])
    else:
        fits_wcs = cube.fits_wcs
        if fits_wcs is None:
            fits_wcs = cube.meta.get("fits_wcs")
    if isinstance(fits_wcs, list):
        fits_wcs = fits_wcs[0] if fits_wcs else None
    if fits_wcs is None or not hasattr(fits_wcs, "celestial"):
        msg = "This cube does not carry the WCS metadata needed to derive a celestial frame."
        raise ValueError(msg)
    return wcs_to_celestial_frame(fits_wcs.celestial)


class _ResolveNegativeIndicesMixin:
    """
    Resolve negative indices before NDCube slicing.

    The sliced WCS keeps negative offsets as they are (astropy#15557), so for example
    ``cube[1:5][-1]`` describes ``cube[0]``, and NDMeta raises on negative slices of
    axis-aware keys.

    TODO: delete this once irispy requires an ndcube that resolves negative indices
    itself (the ndcube ``negative-indices`` branch, which also restores ``meta`` when
    slicing raises).
    """

    def __getitem__(self, item):
        if any(np.iterable(index) for index in (item if isinstance(item, tuple | list) else (item,))):
            # As ndcube, where sanitize_slices rejects an array index, which simplify_basic_index trips over
            msg = "This slice is invalid, only integer or range slices are supported."
            raise IndexError(msg)
        # Against the shape, as since astropy 8.0.2 sanitize_slices rejects negative indices without it
        try:
            item = simplify_basic_index(item, shape=self.shape)
        except np.exceptions.AxisError:
            # An IndexError, as for arrays, not numpy's AxisError, which is a ValueError and speaks of axes
            msg = f"Index {item} is out of range for a cube of shape {self.shape}"
            raise IndexError(msg) from None
        # sanitize_slices raises for any step other than 1, so dropping the step loses nothing.
        sanitize_slices(item, len(self.shape))
        item = tuple(slice(index.start, index.stop) if isinstance(index, slice) else index for index in item)
        if 0 not in self.shape and any(isinstance(index, slice) and index.start >= index.stop for index in item):
            msg = "Slicing would give a length-0 axis, which a WCS cannot describe."
            raise IndexError(msg)
        meta = self.meta
        try:
            return super().__getitem__(item)
        finally:
            # NDCube leaves the unsliced cube without its meta when slicing raises.
            self.meta = meta
