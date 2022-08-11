"""
================
**irispy-lmsal**
================

**irispy-lmsal** is a library that provides the tools to read in and analyze data from Interface Region Imaging Spectrograph (IRIS).
"""

from .sji import *  # NOQA
from .spectrograph import *  # NOQA
from .version import version as __version__

__all__ = ["__version__"]
