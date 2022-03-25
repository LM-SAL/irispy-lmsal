import textwrap

import matplotlib.pyplot as plt
import sunpy.visualization.colormaps as cm  # NOQA
from astropy.time import Time

from sunraster import SpectrogramCube, SpectrogramSequence

from irispy import utils

__all__ = ["IRISMapCube", "IRISMapCubeSequence"]


class IRISMapCube(SpectrogramCube):
    """
    Class representing SJI Image described by a single WCS.

    Parameters
    ----------
    data: `numpy.ndarray`
        The array holding the actual data in this object.
    wcs: `ndcube.wcs.wcs.WCS`
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
        if self.global_coords and "time" in self.global_coords:
            instance_start = self.global_coords["time"].min().isot
            instance_end = self.global_coords["time"].max().isot
        elif self.extra_coords and isinstance(self.axis_world_coords("time", wcs=self.extra_coords)[0], Time):
            instance_start = self.axis_world_coords("time", wcs=self.extra_coords)[0].min().isot
            instance_end = self.axis_world_coords("time", wcs=self.extra_coords)[0].max().isot
        else:
            instance_start = None
            instance_end = None
        return textwrap.dedent(
            f"""
            IRISMapCube
            -----------
            Observatory:\t\t {self.meta.get("TELESCOP")}
            Instrument:\t\t\t {self.meta.get("INSTRUME")}
            Bandpass:\t\t\t {self.meta.get("TWAVE1")}
            Obs. Start:\t\t\t {startobs}
            Obs. End:\t\t\t {endobs}
            Instance Start:\t\t {instance_start}
            Instance End:\t\t {instance_end}
            Roll:\t\t\t {self.meta.get("SAT_ROT")}
            Total Frames in Obs.:\t {self.meta.get("NBFRAMES")}
            IRIS Obs. id:\t\t {self.meta.get("OBSID")}
            IRIS Obs. Description:\t {self.meta.get("OBS_DESC")}
            Axis Types:\t\t\t {self.array_axis_physical_types}
            Cube dimensions:\t\t {self.dimensions}
            """
        )

    def __getitem__(self, item):
        sliced_self = super().__getitem__(item)
        sliced_self.scaled = self.scaled
        return sliced_self

    def plot(self, *args, **kwargs):
        cmap = kwargs.pop("cmap")
        if not cmap:
            cmap = plt.get_cmap(name="irissji{}".format(int(self.meta["TWAVE1"])))
        return super().plot(*args, cmap=cmap, **kwargs)

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
        dust_mask = utils.calculate_dust_mask(self.data)
        if undo:
            # If undo kwarg IS set, unmask dust pixels.
            self.mask[dust_mask] = False
            self.dust_masked = False
        else:
            # If undo kwarg is NOT set, mask dust pixels.
            self.mask[dust_mask] = True
            self.dust_masked = True


class IRISMapCubeSequence(SpectrogramSequence):
    """
    Class for holding, slicing and plotting IRIS SJI data.

    This class contains all the functionality of its super class with
    some additional functionalities.

    Parameters
    ----------
    data_list: `list`
        List of `IRISMapCube` objects from the same OBS ID.
        Each cube must contain the same 'detector type' in its meta attribute.
    meta: `dict` or header object, optional
        Metadata associated with the sequence.
    common_axis: `int`, optional
        The axis of the NDCubes corresponding to time.
    """

    def __init__(self, data_list, meta=None, common_axis=0):
        super().__init__(data_list, meta=meta, common_axis=common_axis)

    def __repr__(self):
        return f"{object.__repr__(self)}\n{str(self)}"

    def __str__(self):
        startobs = self.meta.get("STARTOBS")
        startobs = startobs.isot if startobs else None
        endobs = self.meta.get("ENDOBS")
        endobs = endobs.isot if endobs else None
        instance_start = self[0].extra_coords["time"]["value"]
        instance_start = instance_start.isot if instance_start else None
        instance_end = self[-1].extra_coords["time"]["value"]
        instance_end = instance_end.isot if instance_end else None
        return textwrap.dedent(
            f"""
                IRISMapCubeSequence
                -------------------
                Observatory:\t\t {self.meta.get("TELESCOP")}
                Instrument:\t\t {self.meta.get("INSTRUME")}

                OBS ID:\t\t\t {self.meta.get("OBSID")}
                OBS Description:\t {self.meta.get("OBS_DESC")}
                OBS period:\t\t {startobs} -- {endobs}

                Sequence period:\t {instance_start} -- {instance_end}
                Sequence Shape:\t\t {self.dimensions}
                Roll:\t\t\t {self.meta.get("SAT_ROT")}
                Axis Types:\t\t {self.array_axis_physical_types}
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
