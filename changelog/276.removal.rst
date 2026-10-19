Removed ``irispy.utils.get_detector_type``, which raised `KeyError` on the metadata of every irispy cube.
Use ``cube.meta.detector_band`` (`irispy.meta.BaseMeta.detector_band`) instead, which gives ``"FUV"`` for both FUV detectors.
