"""
This package contains all of irispy's test data.
"""

from pathlib import Path

import irispy

__all__ = [
    "ROOTDIR",
    "get_test_filepath",
]

ROOTDIR = Path(irispy.__file__).parent / "data" / "test"


def get_test_filepath(filename):
    """
    Return the full path to a test file in the ``data/test`` directory.

    Parameters
    ----------
    filename : `str` or `pathlib.Path`
        The name of the file inside the ``data/test`` directory.

    Returns
    -------
    filepath : `str`
        The full path to the file.

    Raises
    ------
    FileNotFoundError
        If the file is not in the ``data/test`` directory.
    """
    filepath = ROOTDIR / filename
    if not filepath.is_file():
        msg = f"No test file {filepath}"
        raise FileNotFoundError(msg)
    return str(filepath)
