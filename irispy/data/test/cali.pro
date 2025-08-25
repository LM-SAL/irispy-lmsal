; Creates the calibration test files
compile_opt idl2
raster_filename = 'iris_l2_20210905_001833_3620258102_raster_t000_r00000_test.fits'
read_iris_l2, raster_filename, index, data
; "C II 1336",
input_spectrum = reform(data[11, 20, *])
wavelength = indgen(187) * 0.025960000231900000 + 1332.6762501600001
save, input_spectrum, wavelength, filename = 'input_calibration.sav'
calib_spectrum = IRIS_CALIB_SPECTRUM(input_spectrum, wavelength, /verbose)
save, calib_spectrum, filename = 'output_calibration.sav'

; This needs the full raster, it does not work with the compressed test version stored in this repo.
raster_filename = 'iris_l2_20210905_001833_3620258102_raster_t000_r00000.fits'
iris_datacal = iris_flux_radcal(local = raster_filename, /onlyfactor, /verbose)
save, iris_datacal, filename = 'iris_datacal.sav'
end
