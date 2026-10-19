`~irispy.obsid.ObsID` now decodes Full readout OBS IDs starting with 41, such as 4180000056, and raises a `ValueError` for OBS IDs that do not match the OBS tables.
For calibration and special sequences, OBS IDs of 4200000000 and above, the error points to the ``OBS_DESC`` and ``OBSTITLE`` header keywords.
