"""
Search and download IRIS data with `sunpy.net.Fido`.

Importing `irispy.net` adds `irispy.net.LMSALClient` to ``Fido`` and the
`irispy.net.attrs` to `sunpy.net.attrs` as ``a.iris``.
"""

from irispy.net.client import LMSALClient

__all__ = ["LMSALClient"]
