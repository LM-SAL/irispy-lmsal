`~irispy.utils.spectrograph.radiometric_calibration` now calibrates slit-jaw cubes (`~irispy.sji.SJICube`).
Its new ``line`` keyword gives the intensity of a single line at that wavelength; the default, `None`, gives a continuum flat across the passband.
