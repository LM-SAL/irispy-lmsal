"""
IRISPy
======

"""
import sys

from .sji import *  # NOQA
from .spectrograph import *  # NOQA

# Enforce Python version check during package import.
__minimum_python_version__ = "3.6"

try:
    from .version import __version__
except ImportError:
    __version__ = "unknown"


class UnsupportedPythonError(Exception):
    pass


if sys.version_info < tuple(int(val) for val in __minimum_python_version__.split(".")):
    # This has to be .format to keep backwards compatibly.
    raise UnsupportedPythonError(
        f"IRISPy does not support Python < {__minimum_python_version__}"
    )
