import pickle

import numpy as np
import pytest

from irispy.io import read_files

# The raster classes also carry a FITS-TAB WCS, which astropy rebuilds from a header alone on unpickling.
NDMETA_ISSUE = "NDMeta cannot be pickled, https://github.com/sunpy/ndcube/issues/847"


@pytest.mark.xfail(strict=True, reason=NDMETA_ISSUE)
@pytest.mark.parametrize(
    ("file_fixture", "keys"),
    [
        ("sns_sji_1330_file", ()),
        ("sot_fg_file", ()),
        ("sns_sg_file", ()),
        ("sns_sg_file", ("Si IV 1403",)),
        ("sns_sg_file", ("Si IV 1403", 0)),
    ],
    ids=["SJICube", "SOTCube", "RasterCollection", "SpectrogramCubeSequence", "SpectrogramCube"],
)
def test_pickle_round_trip(request, file_fixture, keys):
    original = read_files(request.getfixturevalue(file_fixture))
    for key in keys:
        original = original[key]

    copy = pickle.loads(pickle.dumps(original))  # NOQA: S301

    assert type(copy) is type(original)
    if hasattr(original, "data"):
        np.testing.assert_array_equal(copy.data, original.data)
        assert copy.unit == original.unit
    else:
        assert list(copy.keys()) == list(original.keys())
    assert len(copy) == len(original)
