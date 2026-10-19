from pathlib import Path

import irispy

__all__ = ["download_all"]  # NOQA: F822

ROOTDIR = Path(irispy.__file__).parent / "data"


# See PEP 562 (https://peps.python.org/pep-0562/): the sample data downloader is imported on first use
def __getattr__(name):
    if name == "download_all":
        from irispy.data.sample import download_all  # NOQA: PLC0415

        return download_all
    msg = f"module {__name__!r} has no attribute {name!r}"
    raise AttributeError(msg)
