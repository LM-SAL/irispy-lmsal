import os

import irispy
from irispy.data.sample import download_all

__all__ = ["download_all"]

ROOTDIR = os.path.join(os.path.dirname(irispy.__file__), "data")
