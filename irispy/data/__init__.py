from pathlib import Path

import irispy
from irispy.data.sample import download_all

__all__ = ["download_all"]

ROOTDIR = Path(irispy.__file__).parent / "data"
