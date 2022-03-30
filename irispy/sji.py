import textwrap
from collections import defaultdict

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.time import Time
from ndcube.extra_coords import ExtraCoords

from sunraster import SpectrogramCube, SpectrogramSequence

from irispy.utils import calculate_dust_mask
from irispy.visualization import IRISSequencePlotter, _set_axis_colors

__all__ = ["IRISMapCube", "IRISMapCubeSequence"]


def _get_times(iris_map_cube):
    instance_start = None
    instance_end = None
    if hasattr(iris_map_cube, "global_coords") and "time" in iris_map_cube.global_coords:
        instance_start = iris_map_cube.global_coords["time"].min().isot
        instance_end = iris_map_cube.global_coords["time"].max().isot
    elif hasattr(iris_map_cube, "time") and iris_map_cube.time:
        instance_start = iris_map_cube.time.min().isot
        instance_end = iris_map_cube.time.max().isot
    elif (
        hasattr(iris_map_cube, "extra_coords")
        and hasattr(iris_map_cube, "axis_world_coords")
        and isinstance(iris_map_cube.axis_world_coords("time", wcs=iris_map_cube.extra_coords)[0], Time)
    ):
        instance_start = iris_map_cube.axis_world_coords("time", wcs=iris_map_cube.extra_coords)[0].min().isot
        instance_end = iris_map_cube.axis_world_coords("time", wcs=iris_map_cube.extra_coords)[0].max().isot
    return instance_start, instance_end


class IRISMapCube(SpectrogramCube):
    """
    Class representing SJI Image described by a single WCS.

    Parameters
    ----------
    data : `numpy.ndarray`
        The array holding the actual data in this object.
    wcs : `ndcube.wcs.wcs.WCS`
        The WCS object containing the axes' information
    unit : `astropy.unit.Unit` or `str`
        Unit for the dataset.
        Strings that can be converted to a Unit are allowed.
    meta : dict-like object
        Additional meta information about the dataset.
    uncertainty : any type, optional
        Uncertainty in the dataset. Should have an attribute uncertainty_type
        that defines what kind of uncertainty is stored, for example "std"
        for standard deviation or "var" for variance. A metaclass defining
        such an interface is NDUncertainty - but isn't mandatory. If the
        uncertainty has no such attribute the uncertainty is stored as
        UnknownUncertainty.
        Defaults to None.
    mask : any type, optional
        Mask for the dataset. Masks should follow the numpy convention
        that valid data points are marked by False and invalid ones with True.
        Defaults to None.
    copy : `bool`, optional
        Indicates whether to save the arguments as copy. True copies every
        attribute before saving it while False tries to save every parameter
        as reference. Note however that it is not always possible to save the
        input as reference.
        Default is False.
    scaled : `bool`, optional
        Indicates if the data has been scaled.
    """

    def __init__(
        self,
        data,
        wcs,
        uncertainty=None,
        unit=None,
        meta=None,
        mask=None,
        copy=False,
        scaled=None,
        **kwargs,
    ):
        self.scaled = scaled
        self.dust_masked = False
        super().__init__(
            data,
            wcs,
            uncertainty=uncertainty,
            mask=mask,
            meta=meta,
            unit=unit,
            copy=copy,
            **kwargs,
        )

    def __repr__(self):
        return f"{object.__repr__(self)}\n{str(self)}"

    def __str__(self):
        startobs = self.meta.get("STARTOBS")
        startobs = startobs.isot if startobs else None
        endobs = self.meta.get("ENDOBS")
        endobs = endobs.isot if endobs else None
        instance_start, instance_end = _get_times(self)
        return textwrap.dedent(
            f"""
            IRISMapCube
            -----------
            Observatory:           {self.meta.get("TELESCOP")}
            Instrument:            {self.meta.get("INSTRUME")}
            Bandpass:              {self.meta.get("TWAVE1")}
            Obs. Start:            {startobs}
            Obs. End:              {endobs}
            Instance Start:        {instance_start}
            Instance End:          {instance_end}
            Total Frames in Obs.:  {self.meta.get("NBFRAMES")}
            IRIS Obs. id:          {self.meta.get("OBSID")}
            IRIS Obs. Description: {self.meta.get("OBS_DESC")}
            Axis Types:            {self.array_axis_physical_types}
            Roll:                  {self.meta.get("SAT_ROT")}
            Cube dimensions:       {self.dimensions}
            """
        )

    def __getitem__(self, item):
        sliced_self = super().__getitem__(item)
        sliced_self.scaled = self.scaled
        return sliced_self

    def plot(self, *args, **kwargs):
        cmap = kwargs.get("cmap")
        if not cmap:
            try:
                cmap = plt.get_cmap(name="irissji{}".format(int(self.meta["TWAVE1"])))
            except Exception:
                cmap = "viridis"
        kwargs["cmap"] = cmap
        ax = super().plot(*args, **kwargs)
        _set_axis_colors(ax)
        return ax

    def apply_dust_mask(self, undo=False):
        """
        Applies or undoes an update of the mask with the dust particles
        positions.

        Rewrite self.mask with/without the dust positions.

        Parameters
        ----------
        undo: `bool`
            If False, dust particles positions mask will be applied.
            If True, dust particles positions mask will be removed.
            Default=False
        """
        dust_mask = calculate_dust_mask(self.data)
        if undo:
            # If undo kwarg IS set, unmask dust pixels.
            self.mask[dust_mask] = False
            self.dust_masked = False
        else:
            # If undo kwarg is NOT set, mask dust pixels.
            self.mask[dust_mask] = True
            self.dust_masked = True


def _create_extra_coords(cube):
    extra_coords = ExtraCoords()
    # = [("time", 0, times)]
    coords = defaultdict(list)
    for data in cube.data:
        for name, value in data.global_coords.items():
            coords[name].append(value)
    for k, v in coords.items():
        if k == "time":
            v = Time(v)
        else:
            v = u.Quantity(v)
        extra_coords.add(k, 0, v)
    # extra_coords.add([name, value for name, value in coords.items()])
    return extra_coords


class IRISMapCubeSequence(SpectrogramSequence):
    """
    Class for holding, slicing and plotting IRIS SJI data.

    This class contains all the functionality of its super class with
    some additional functionalities.

    Parameters
    ----------
    data_list: `list`
        List of `IRISMapCube` objects from the same OBS ID.
    meta: `dict` or header object, optional
        Metadata associated with the sequence.
    common_axis: `int`, optional
        The axis of the NDCubes corresponding to time.
    """

    plotter = IRISSequencePlotter

    def __init__(self, data_list, meta=None, common_axis=0, times=None):
        super().__init__(data_list, meta=meta, common_axis=common_axis)
        self.time = times
        self._extra_coords = _create_extra_coords(self)
        self._data = np.stack(list(data.data for data in data_list))

    def __repr__(self):
        return f"{object.__repr__(self)}\n{str(self)}"

    def __str__(self):
        startobs = self.meta.get("STARTOBS")
        startobs = startobs.isot if startobs else None
        endobs = self.meta.get("ENDOBS")
        endobs = endobs.isot if endobs else None
        instance_start, instance_end = _get_times(self)
        return textwrap.dedent(
            f"""
            IRISMapCubeSequence
            -------------------
            Observatory:     {self.meta.get("TELESCOP")}
            Instrument:      {self.meta.get("INSTRUME")}
            OBS ID:          {self.meta.get("OBSID")}
            OBS Description: {self.meta.get("OBS_DESC")}
            OBS period:      {startobs} -- {endobs}
            Sequence period: {instance_start} -- {instance_end}
            Sequence Shape:  {self.dimensions}
            Axis Types:      {self.array_axis_physical_types}
            Roll:            {self.meta.get("SAT_ROT")}
            """
        )

    def apply_dust_mask(self, undo=False):
        """
        Applies or undoes an update of all the masks with the dust particles
        positions.

        Rewrites all data with/without the dust positions.

        Parameters
        ----------
        undo: `bool`, optional
            If `False` (default), dust particles positions masks will be applied.
            If `True`, dust particles positions masks will be removed.
        """
        for cube in self.data:
            cube.apply_dust_mask(undo=undo)

    @property
    def time(self):
        return self._time

    @time.setter
    def time(self, times):
        self._time = times

    @property
    def exposure_time(self):
        return self._extra_coords["exposure time"]

    @property
    def extra_coords(self):
        return self._extra_coords

    @property
    def data_as_array(self):
        return self._data

    def plot(self, *args, **kwargs):
        cmap = kwargs.get("cmap")
        if not cmap:
            try:
                cmap = plt.get_cmap(name="irissji{}".format(int(self.meta["TWAVE1"])))
            except Exception:
                cmap = "viridis"
        kwargs["cmap"] = cmap
        ax = super().plot(self, *args, **kwargs)
        _set_axis_colors(ax)
        return ax
