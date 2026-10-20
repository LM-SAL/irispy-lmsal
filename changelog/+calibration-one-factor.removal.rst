Removed ``irispy.utils.spectrograph.convert_photons_per_sec_to_radiance``; multiply data in photons per second by ``irispy.utils.spectrograph.calculate_dn_to_radiance_factor`` instead.
Removed ``irispy.utils.spectrograph.reshape_1d_wavelength_dimensions_for_broadcast``; an array of wavelengths already broadcasts along the last axis of the data.
