``irispy.obsid.ObsID`` no longer leaves trailing spaces in the ``options["sjis"]`` keys, so ``decoded.options["sjis"][decoded["sjis"]]`` works instead of raising a ``KeyError``.
