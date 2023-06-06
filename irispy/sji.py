import logging
import textwrap

import matplotlib.pyplot as plt
from sunraster import SpectrogramCube

from irispy.utils import calculate_dust_mask
from irispy.visualization import Plotter, _set_axis_colors

__all__ = ["SJICube"]


class SJICube(SpectrogramCube):
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
        *,
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
        if "_basic_wcs" in kwargs:
            _basic_wcs = kwargs.pop("_basic_wcs")
        else:
            _basic_wcs = None
        self._basic_wcs = _basic_wcs
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
        if self.wcs.world_n_dim == 2:
            instance_start = self.global_coords["Time (UTC)"]
            instance_end = None
        else:
            instance_start = self.wcs.pixel_to_world(0, 0, 0)[-1]
            instance_end = self.wcs.pixel_to_world(0, 0, self.data.shape[0] - 1)[-1]
        return textwrap.dedent(
            f"""
            SJICube
            -------
            Observatory:           {self.meta.get("TELESCOP")}
            Instrument:            {self.meta.get("INSTRUME")}
            Bandpass:              {self.meta.get("TWAVE1")}
            Obs. Start:            {self.meta.get("STARTOBS")}
            Obs. End:              {self.meta.get("ENDOBS")}
            Instance Start:        {instance_start}
            Instance End:          {instance_end}
            Total Frames in Obs.:  {self.meta.get("NBFRAMES")}
            IRIS Obs. id:          {self.meta.get("OBSID")}
            IRIS Obs. Description: {self.meta.get("OBS_DESC")}
            Axis Types:            {self.array_axis_physical_types}
            Roll:                  {self.meta.get("SAT_ROT")}
            Cube dimensions:       {self.dimensions}
            """,
        )

    def __getitem__(self, item):
        sliced_self = super().__getitem__(item)
        sliced_self.scaled = self.scaled
        if self._basic_wcs is not None:
            sliced_self._basic_wcs = self._basic_wcs[item]
        return sliced_self

    def plot(self, *args, **kwargs):
        cmap = kwargs.get("cmap")
        if not cmap:
            try:
                cmap = plt.get_cmap(name=f"irissji{int(self.meta['TWAVE1'])}")
            except Exception as e:  # NOQA: BLE001
                logging.debug(e)
                cmap = "viridis"
        kwargs["cmap"] = cmap
        ax = Plotter(ndcube=self).plot(*args, **kwargs)
        _set_axis_colors(ax)
        return ax

    def apply_dust_mask(self, *, undo=False):
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

    @property
    def basic_wcs(self):
        """
        Returns a standard WCS instead of gWCS.
        """
        return self._basic_wcs
