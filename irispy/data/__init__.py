import os

import irispy
from irispy.data._sample import download_sample_data

__all__ = ["download_sample_data"]

rootdir = os.path.join(os.path.dirname(irispy.__file__), "data")
