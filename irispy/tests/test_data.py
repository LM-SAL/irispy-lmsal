import pytest

from irispy.data.test import get_test_filepath


def test_get_test_filepath_missing_file():
    with pytest.raises(FileNotFoundError, match=r"iris_l2_typo\.fits"):
        get_test_filepath("iris_l2_typo.fits")
