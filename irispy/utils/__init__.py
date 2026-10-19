"""
Analysis and calibration helpers.

The submodules are imported on first use, so ``import irispy`` does not pay for them.
"""

from importlib import import_module

from .utils import *  # This one is to avoid a double utils in import path
from .utils import __all__ as _utils_all

_SUBMODULES = (
    "bursts",
    "constants",
    "cosmic_rays",
    "density",
    "dust",
    "fiducials",
    "fitting",
    "mg_features",
    "moments",
    "red_blue",
    "response",
    "rgb",
    "spectrograph",
    "wavelength_drift",
    "wobble",
)
__all__ = [*_SUBMODULES, *_utils_all, "calculate_dust_mask"]  # NOQA: PLE0604, F405


# See PEP 562 (https://peps.python.org/pep-0562/) for module-level __getattr__() and __dir__()
def __getattr__(name):
    if name in _SUBMODULES:
        return import_module(f"{__name__}.{name}")
    if name == "calculate_dust_mask":
        from .dust import calculate_dust_mask  # NOQA: PLC0415

        return calculate_dust_mask
    msg = f"module {__name__!r} has no attribute {name!r}"
    raise AttributeError(msg)


def __dir__():
    # utils is the submodule behind the star import above
    return sorted([*__all__, "utils"])
