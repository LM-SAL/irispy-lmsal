"""
Utilities for repairing dust-darkened pixels in IRIS image cubes.
"""

import warnings
import functools

import numpy as np
from scipy import ndimage

import astropy.units as u
from astropy.io import fits

from irispy.data import ROOTDIR
from ._spectral import check_scaled
from .utils import _boolean_mask

__all__ = ["calculate_dust_mask", "remove_dust"]


@functools.cache
def _badpix_table():
    """
    The detector pixels and pointing constants per channel, and the level 2 plate scale,
    from the newest pipeline bad-pixel table in ``irispy/data``.
    """
    path = max(ROOTDIR.glob("iris_sji_badpix_*.fits"))
    with fits.open(path) as hdul:
        pointing = {
            str(row["CHANNEL"]).strip(): (row["CPX1"], row["CPX2"], row["CDELT"], row["ROLL"])
            for row in hdul["POINTING"].data
        }
        pixels = {
            channel: (hdul[channel].data["X"].astype(int), hdul[channel].data["Y"].astype(int)) for channel in pointing
        }
        return pixels, pointing, hdul[0].header["L2CDELT"]


def calculate_dust_mask(cube):
    """
    The IRIS pipeline's bad pixels (dust and dead pixels) in the frame of a slit-jaw
    cube.

    Parameters
    ----------
    cube : `~irispy.sji.SJICube`
        A slit-jaw cube read from a level 2 file, or a slice of one in time or space.

    Returns
    -------
    `numpy.ndarray` of `bool`
        `True` at the bad pixels and their neighbors, with the shape of the cube's
        data. Pixels without a value (not finite) stay `False`.
    """
    channel = cube.meta.get("TDESC1")
    pixels, pointing, level2_cdelt = _badpix_table()
    if channel not in pixels:
        msg = f"No bad-pixel table for {channel!r}; the tables cover {sorted(pixels)}."
        raise ValueError(msg)
    if cube._file_frames is None:
        msg = "The bad pixels can be placed only in a cube read from a level 2 slit-jaw file and not rebinned."
        raise ValueError(msg)
    spatial = cube._file_frames[1]
    if not all(isinstance(axis, slice) for axis in spatial):
        msg = "The bad pixels can be placed only in an image, not in a 1-D cut of one."
        raise ValueError(msg)
    # IRIS_DUSTBUSTER (SolarSoft): a bad pixel's summed detector pixel relative to the slit center
    # (the pointing database's 1-based FITS pixel, at the center of its summed pixel), rotated by
    # the detector's roll against the slit and scaled from the detector's plate scale to the
    # level 2 one; -0.5 is the routine's default y adjustment.
    detector_x, detector_y = pixels[channel]
    cpx1, cpx2, cdelt, roll = pointing[channel]
    sum_x, sum_y = cube.meta["SUMSPTRL"], cube.meta["SUMSPAT"]
    dx = detector_x // sum_x - ((cpx1 - 1) / sum_x + (sum_x - 1) / (2 * sum_x))
    dy = detector_y // sum_y - 0.5 - ((cpx2 - 1) / sum_y + (sum_y - 1) / (2 * sum_y))
    distance = np.hypot(dx, dy) * cdelt / level2_cdelt
    angle = np.arctan2(dx, dy) - np.deg2rad(roll)
    offset_x = distance * np.sin(angle) - (spatial[1].start or 0)
    offset_y = distance * np.cos(angle) - (spatial[0].start or 0)
    if cube.data.ndim == 2:
        slit_x = np.atleast_1d(cube.global_coords["slit x position"].to_value(u.pix))
        slit_y = np.atleast_1d(cube.global_coords["slit y position"].to_value(u.pix))
    else:
        slit_x = cube.axis_world_coords(wcs=cube.extra_coords["slit x position"])[0].to_value(u.pix)
        slit_y = cube.axis_world_coords(wcs=cube.extra_coords["slit y position"])[0].to_value(u.pix)
    ny, nx = cube.data.shape[-2:]
    mask = np.zeros((len(slit_x), ny, nx), dtype=bool)
    x = offset_x + slit_x[:, np.newaxis]
    y = offset_y + slit_y[:, np.newaxis]
    frame = np.broadcast_to(np.arange(len(slit_x))[:, np.newaxis], x.shape)
    # Floor and ceil in both axes, as the IDL routine does, so the mapped mask has no gaps
    for xi in (np.floor(x), np.ceil(x)):
        for yi in (np.floor(y), np.ceil(y)):
            inside = (xi >= 0) & (xi < nx) & (yi >= 0) & (yi < ny)
            mask[frame[inside], yi[inside].astype(int), xi[inside].astype(int)] = True
    # The list is one level 2 pixel tight on the images; IRIS_DUSTBUSTER smooths it over 5 detector pixels
    mask = ndimage.binary_dilation(mask, structure=np.ones((1, 3, 3), dtype=bool)).reshape(cube.data.shape)
    # The dead columns of a record run through the frame's edge, which holds no values
    return mask & np.asarray(np.isfinite(cube.data))


def _coerce_dust_mask(cube, data, dust_mask):
    dust_mask = calculate_dust_mask(cube) if dust_mask is None else np.asarray(dust_mask, dtype=bool)
    if data.ndim == 3 and dust_mask.shape == data.shape[1:]:
        dust_mask = np.broadcast_to(dust_mask, data.shape)
    if dust_mask.shape != data.shape:
        msg = f"dust_mask must have shape {data.shape} or {data.shape[1:]}; got {dust_mask.shape}"
        raise ValueError(msg)
    # Pixels without a value are not dust: filling them would paint over the frame's edges
    return dust_mask & np.isfinite(data)


def _local_median_fill(frame, invalid_mask, spatial_box, target_mask):
    """
    Local ``nanmedian`` of ``frame`` at ``target_mask`` pixels, NaN elsewhere.

    Matches ``generic_filter`` with ``np.nanmedian`` and ``mode="nearest"`` at the
    target pixels, without filtering the whole frame.
    """
    masked_frame = np.where(invalid_mask | ~np.isfinite(frame), np.nan, frame)
    padded = np.pad(masked_frame, spatial_box // 2, mode="edge")
    windows = np.lib.stride_tricks.sliding_window_view(padded, (spatial_box, spatial_box))
    rows, cols = np.nonzero(target_mask)
    spatial_fill = np.full(frame.shape, np.nan)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        spatial_fill[rows, cols] = np.nanmedian(windows[rows, cols].reshape(rows.size, spatial_box**2), axis=1)
    return spatial_fill


def _resolve_exposure_times(cube, frame_count, *, exposure_normalize):
    """
    Resolve per-frame exposure times for temporal normalization.

    Returns
    -------
    `numpy.ndarray` or None
        Per-frame exposure times in seconds when available and valid.
        Returns None when normalization is disabled or metadata are unusable.
    """
    if not exposure_normalize:
        return None
    try:
        exposure_times = u.Quantity(cube.exposure_time, copy=False).to_value(u.s)
    except (AttributeError, TypeError, ValueError):
        return None
    exposure_times = np.atleast_1d(np.asarray(exposure_times, dtype=float))
    if exposure_times.size == 1:
        return np.repeat(exposure_times, frame_count)
    if exposure_times.size != frame_count:
        warnings.warn(
            "exposure_normalize=True but the number of exposure_time values "
            f"({exposure_times.size}) does not match the number of frames "
            f"in the cube ({frame_count}). Temporal exposure "
            "normalization has been skipped.",
            UserWarning,
            stacklevel=3,
        )
        return None
    return exposure_times


def remove_dust(
    cube,
    *,
    dust_mask=None,
    temporal_window=2,
    exposure_normalize=True,
    fallback="spatial",
    spatial_box=5,
):
    """
    Replace dust-darkened pixels in an IRIS image cube.

    This routine is inspired by SolarSoft's ``IRIS_DUSTBUSTER``. It uses
    neighboring frames at the same pixel location first, and falls back to a local
    spatial median when temporal replacements is not available.

    Parameters
    ----------
    cube : `~irispy.sji.SJICube`
        The image cube to clean. Two-dimensional image slices are also supported.
    dust_mask : `numpy.ndarray`, optional
        Boolean mask marking the pixels to repair. If omitted, the pipeline's bad
        pixels from `~irispy.utils.dust.calculate_dust_mask` are repaired.
        For a 3D cube, a 2D mask is broadcast over time.
    temporal_window : `int`, optional
        Number of neighboring frames on either side to consider when computing the
        replacement value. Defaults to 2, which mirrors the neighboring frames used
        in the IDL implementation.
    exposure_normalize : `bool`, optional
        If `True`, normalize candidate replacement pixels by the frame exposure time
        before taking the median and scale the result back to the target frame.
        If no usable exposure time metadata are available this step is skipped.
    fallback : {``"spatial"``, None}, optional
        Fallback to use when temporal replacement is unavailable.
        ``"spatial"`` applies a local median filter to the frame.
        None leaves those pixels masked.
    spatial_box : `int`, optional
        Size of the local median filter used by the spatial fallback.
        Defaults to 5.

    Returns
    -------
    `~irispy.sji.SJICube`
        A new cube with repaired dust pixels.

    Notes
    -----
    Pixels without a value (not finite) are never repaired.
    """
    if fallback not in {"spatial", None}:
        msg = "fallback must be 'spatial' or None."
        raise ValueError(msg)

    check_scaled(cube)
    clean_data = np.asarray(cube.data, dtype=float).copy()
    original_mask = _boolean_mask(cube)
    dust_mask = _coerce_dust_mask(cube, clean_data, dust_mask)
    filled_mask = np.zeros_like(dust_mask, dtype=bool)

    if clean_data.ndim == 3 and temporal_window > 0 and np.any(dust_mask):
        exposure_times = _resolve_exposure_times(
            cube,
            clean_data.shape[0],
            exposure_normalize=exposure_normalize,
        )
        for frame_index in range(clean_data.shape[0]):
            target_mask = dust_mask[frame_index]
            if not np.any(target_mask):
                continue
            candidate_frames = []
            for offset in range(1, temporal_window + 1):
                for neighbor_index in (frame_index - offset, frame_index + offset):
                    if 0 <= neighbor_index < clean_data.shape[0]:
                        candidate = clean_data[neighbor_index].copy()
                        invalid = original_mask[neighbor_index] | dust_mask[neighbor_index] | ~np.isfinite(candidate)
                        candidate[invalid] = np.nan
                        if exposure_times is not None:
                            candidate = candidate / exposure_times[neighbor_index] * exposure_times[frame_index]
                        candidate_frames.append(candidate)
            if candidate_frames:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", category=RuntimeWarning)
                    temporal_fill = np.nanmedian(np.stack(candidate_frames), axis=0)
                can_fill = target_mask & np.isfinite(temporal_fill)
                clean_data[frame_index, can_fill] = temporal_fill[can_fill]
                filled_mask[frame_index, can_fill] = True

    if fallback == "spatial":
        if clean_data.ndim == 2:
            spatial_fill = _local_median_fill(clean_data, original_mask | dust_mask, spatial_box, dust_mask)
            can_fill = dust_mask & np.isfinite(spatial_fill)
            clean_data[can_fill] = spatial_fill[can_fill]
            filled_mask[can_fill] = True
        else:
            for frame_index in range(clean_data.shape[0]):
                target_mask = dust_mask[frame_index] & ~filled_mask[frame_index]
                if not np.any(target_mask):
                    continue
                spatial_fill = _local_median_fill(
                    clean_data[frame_index],
                    original_mask[frame_index] | dust_mask[frame_index],
                    spatial_box,
                    target_mask,
                )
                can_fill = target_mask & np.isfinite(spatial_fill)
                clean_data[frame_index, can_fill] = spatial_fill[can_fill]
                filled_mask[frame_index, can_fill] = True

    output_mask = original_mask.copy()
    output_mask[dust_mask] = True
    output_mask[filled_mask] = False
    cleaned_cube_kwargs = {
        "data": clean_data,
        "mask": output_mask,
        "nddata_type": type(cube),
        "extra_coords": "copy",
        "global_coords": "copy",
    }
    cleaned_cube = cube.to_nddata(**cleaned_cube_kwargs)
    if hasattr(cleaned_cube, "dust_masked"):
        cleaned_cube.dust_masked = bool(np.any(output_mask & dust_mask))
    return cleaned_cube
