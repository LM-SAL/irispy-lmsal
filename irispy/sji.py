"""
This module provides movie tools for level 2 IRIS SJI fits file.
"""
import textwrap

import matplotlib.pyplot as plt
import sunpy.visualization.colormaps as cm  # NOQA
from astropy.time import Time

from sunraster import SpectrogramCube, SpectrogramSequence

from irispy import utils

__all__ = ["IRISMapCube", "IRISMapCubeSequence"]


class IRISMapCube(SpectrogramCube):
    """
    Class representing SJI Images described by a single WCS.

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

    Examples
    --------
    >>> from irispy.io import read_sji_lvl2  # doctest: +SKIP
    >>> from irispy.data import sample  # doctest: +SKIP
    >>> sji = read_iris_sji_level2_fits(sample.SJI_CUBE_1400)  # doctest: +SKIP
    >>> sji  # doctest: +SKIP
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
        """
        Initialization of Slit Jaw Imager.
        """
        # Set whether SJI data is scaled or not.
        self.scaled = scaled
        # Dust_masked variable shows whether the dust pixels are set to True in the data mask.
        self.dust_masked = False
        super().__init__(
            data,
            wcs,
            uncertainty=uncertainty,
            mask=mask,
            meta=meta,
            unit=unit,
            copy=copy,
        )

    def __repr__(self):
        return f"{object.__repr__(self)}\n{str(self)}"

    def __str__(self):
        roll = self.meta.get("SAT_ROT", None)
        # Conversion of the start date of OBS
        startobs = self.meta.get("STARTOBS", None)
        startobs = startobs.isot if startobs else None
        # Conversion of the end date of OBS
        endobs = self.meta.get("ENDOBS", None)
        endobs = endobs.isot if endobs else None
        # Conversion of the instance start and end of OBS
        if self.global_coords and "time" in self.global_coords:
            instance_start = self.global_coords["time"].min().isot
            instance_end = self.global_coords["time"].max().isot
        elif self.extra_coords and isinstance(self.axis_world_coords("time", wcs=self.extra_coords)[0], Time):
            instance_start = self.axis_world_coords("time", wcs=self.extra_coords)[0].min().isot
            instance_end = self.axis_world_coords("time", wcs=self.extra_coords)[0].max().isot
        else:
            instance_start = None
            instance_end = None
        # Representation of IRISMapCube object
        return textwrap.dedent(
            f"""
            IRISMapCube
            -----------
            Observatory:\t\t {self.meta.get("TELESCOP", None)}
            Instrument:\t\t\t {self.meta.get("INSTRUME", None)}
            Bandpass:\t\t\t {self.meta.get("TWAVE1", None)}
            Obs. Start:\t\t\t {startobs}
            Obs. End:\t\t\t {endobs}
            Instance Start:\t\t {instance_start}
            Instance End:\t\t {instance_end}
            Roll:\t\t\t {roll}
            Total Frames in Obs.:\t {self.meta.get("NBFRAMES", None)}
            IRIS Obs. id:\t\t {self.meta.get("OBSID", None)}
            IRIS Obs. Description:\t {self.meta.get("OBS_DESC", None)}
            Axis Types:\t\t\t {self.array_axis_physical_types}
            Cube dimensions:\t\t {self.dimensions}
            """
        )

    def __getitem__(self, item):
        # Use parent method to slice item.
        sliced_self = super().__getitem__(item)
        # Ensure value of scaled property is maintained.
        sliced_self.scaled = self.scaled
        return sliced_self

    def plot(self, *args, **kwargs):
        # If colormap not set, load one default sunpy colormap based on SJI passband.
        cmap = kwargs.pop("cmap", None)
        if not cmap:
            cmap = plt.get_cmap(name="irissji{}".format(int(self.meta["TWAVE1"])))

        # Call parent plot function
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
        # Calculate position of dust pixels
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
        # Initialize IRISMapCubeSequence.
        super().__init__(data_list, meta=meta, common_axis=common_axis)

    def __repr__(self):
        return f"{object.__repr__(self)}\n{str(self)}"

    def __str__(self):
        roll = self.meta.get("SAT_ROT", None)
        # Conversion of the start date of OBS
        startobs = self.meta.get("STARTOBS", None)
        startobs = startobs.isot if startobs else None
        # Conversion of the end date of OBS
        endobs = self.meta.get("ENDOBS", None)
        endobs = endobs.isot if endobs else None
        # Conversion of the instance start of OBS
        instance_start = self[0].extra_coords["time"]["value"]
        instance_start = instance_start.isot if instance_start else None
        # Conversion of the instance end of OBS
        instance_end = self[-1].extra_coords["time"]["value"]
        instance_end = instance_end.isot if instance_end else None
        # Representation of IRISMapCube object
        return textwrap.dedent(
            f"""
                IRISMapCubeSequence
                -------------------
                Observatory:\t\t {self.meta.get("TELESCOP", None)}
                Instrument:\t\t {self.meta.get("INSTRUME", None)}

                OBS ID:\t\t\t {self.meta.get("OBSID", None)}
                OBS Description:\t {self.meta.get("OBS_DESC", None)}
                OBS period:\t\t {startobs} -- {endobs}

                Sequence period:\t {instance_start} -- {inst_end}
                Sequence Shape:\t\t {seq_shape}
                Roll:\t\t\t {roll}
                Axis Types:\t\t {axis_types}

                """.format(
                obs=self.meta.get("TELESCOP", None),
                instrument=self.meta.get("INSTRUME", None),
                obs_id=self.meta.get("OBSID", None),
                obs_desc=self.meta.get("OBS_DESC", None),
                obs_start=startobs,
                obs_end=endobs,
                inst_start=instance_start,
                inst_end=instance_end,
                seq_shape=self.dimensions,
                roll=roll,
                axis_types=self.array_axis_physical_types,
            )
        )

    def __getitem__(self, item):
        return self.index_as_cube[item]

    @property
    def dimensions(self):
        return self.cube_like_dimensions

    def plot(self, axes=None, plot_axis_indices=None, axes_coordinates=None, axes_units=None, data_unit=None, **kwargs):
        """
        Visualizes data in the IRISMapCubeSequence with the sequence axis
        folded into the common axis.

        Based on the cube-like dimensionality of the sequence and value of plot_axis_indices
        kwarg, a Line/Image Plot/Animation is produced.

        Parameters
        ----------
        axes: `astropy.visualization.wcsaxes.core.WCSAxes` or ??? or None.
            The axes to plot onto. If None the current axes will be used.
        plot_axis_indices: `int` or iterable of one or two `int`.
            If two axis indices are given, the sequence is visualized as an image or
            2D animation, assuming the sequence has at least 2 cube-like dimensions.
            The cube-like dimension indicated by the 0th element of plot_axis indices is
            displayed on the x-axis while the cube-like dimension indicated by the 1st
            element of plot_axis_indices is displayed on the y-axis. If only one axis
            index is given (either as an int or a list of one int), then a 1D line
            animation is produced with the indicated cube-like dimension on the x-axis
            and other cube-like dimensions represented by animations sliders.
            Default=[-1, -2]. If sequence only has one cube-like dimension,
            plot_axis_indices is ignored and a static 1D line plot is produced.
        axes_coordinates: `None` or `list` of `None` `astropy.units.Quantity` `numpy.ndarray` `str`
            Denotes physical coordinates for plot and slider axes.
            If None coordinates derived from the WCS objects will be used for all axes.
            If a list, its length should equal either the number cube-like dimensions or
            the length of plot_axis_indices.
            If the length equals the number of cube-like dimensions, each element describes
            the coordinates of the corresponding cube-like dimension.
            If the length equals the length of plot_axis_indices,
            the 0th entry describes the coordinates of the x-axis
            while (if length is 2) the 1st entry describes the coordinates of the y-axis.
            Slider axes are implicitly set to None.
            If the number of cube-like dimensions equals the length of plot_axis_indices,
            the latter convention takes precedence.
            The value of each entry should be either
            `None` (implies derive the coordinates from the WCS objects),
            an `astropy.units.Quantity` or a `numpy.ndarray` of coordinates for each pixel,
            or a `str` denoting a valid extra coordinate.
        axes_units: `None or `list` of `None`, `astropy.units.Unit` and/or `str`
            If None units derived from the WCS objects will be used for all axes.
            If a list, its length should equal either the number cube-like dimensions or
            the length of plot_axis_indices.
            If the length equals the number of cube-like dimensions, each element gives the
            unit in which the coordinates along the corresponding cube-like dimension should
            displayed whether they be a plot axes or a slider axes.
            If the length equals the length of plot_axis_indices,
            the 0th entry describes the unit in which the x-axis coordinates should be displayed
            while (if length is 2) the 1st entry describes the unit in which the y-axis should
            be displayed.  Slider axes are implicitly set to None.
            If the number of cube-like dimensions equals the length of plot_axis_indices,
            the latter convention takes precedence.
            The value of each entry should be either
            `None` (implies derive the unit from the WCS object of the 0th sub-cube),
            `astropy.units.Unit` or a valid unit `str`.
        data_unit: `astropy.unit.Unit` or valid unit `str` or None
            Unit in which data be displayed.  If the length of plot_axis_indices is 2,
            a 2D image/animation is produced and data_unit determines the unit represented by
            the color table.  If the length of plot_axis_indices is 1,
            a 1D plot/animation is produced and data_unit determines the unit in which the
            y-axis is displayed.

        Returns
        -------
        `ndcube.mixins.sequence_plotting.plot_as_cube`
        """
        return self.plot_as_cube(
            axes=axes,
            plot_axis_indices=plot_axis_indices,
            axes_coordinates=axes_coordinates,
            axes_units=axes_units,
            data_unit=data_unit,
            **kwargs,
        )

    def apply_dust_mask(self, undo=False):
        """
        Applies or undoes an update of all the masks with the dust particles
        positions.

        Rewrites all self.data[i] mask with/without the dust positions.

        Parameters
        ----------
        undo: `bool`, optional
            If False, dust particles positions masks will be applied.
            If True, dust particles positions masks will be removed.
            Default=False
        """
        for cube in self.data:
            cube.apply_dust_mask(undo=undo)
