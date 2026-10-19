; Creates the IDL reference tables for the non-thermal velocity tests in irispy/utils/tests/test_fitting.py
; (see overview.txt). It reads no data. Needs SolarSoft and this path, as iris_nonthermalwidth reads
; eis_width2velocity.dat and calls eis_element2mass from Hinode/EIS:
;   !path = expand_path('+~/ssw/iris/idl') + ':' + expand_path('+~/ssw/gen/idl') + ':' $
;     + expand_path('~/ssw/hinode/eis/idl/atest/hwarren') + ':' + !path

; The shortest string that reads back as the same double, which is how astropy writes floats.
function shortest, value
  compile_opt idl2
  for n = 1, 17 do begin
    s = string(double(value), format='(F0.' + strtrim(n, 2) + ')')
    if double(s) eq double(value) then return, s
  endfor
end

pro make_idl_references
  compile_opt idl2
  ; iris_nonthermalwidth's speed of light, to give it the observed 1/e width in km/s.
  cc = 2.9979d5
  elements = ['Si', 'O', 'C', 'Fe', 'Fe', 'Mg']
  ions = ['IV', 'IV', 'II', 'XII', 'XXI', 'II']
  wavelengths = [1402.77d, 1401.16d, 1335.71d, 1349.40d, 1354.08d, 2796.35d]
  ; irispy's INSTRUMENTAL_FWHM and the line catalog's log_t_max.
  instrumental = [0.026d, 0.026d, 0.026d, 0.026d, 0.026d, 0.053d]
  log_t = [4.9d, 5.15d, 4.4d, 6.15d, 7.05d, 4.15d]
  ; The first width of each line is below the instrumental and thermal widths added in quadrature.
  fwhm = [[0.05d, 0.08d, 0.15d, 0.3d], [0.08d, 0.12d, 0.2d, 0.3d], [0.04d, 0.08d, 0.15d, 0.3d], $
    [0.12d, 0.2d, 0.3d, 0.5d], [0.4d, 0.6d, 0.9d, 1.2d], [0.06d, 0.15d, 0.3d, 0.5d]]
  ; The instrumental FWHM in the header of iris_nonthermalwidth: 31.8 mA for FUV long (FUV2), 28.6 mA for FUV short
  ; (FUV1). It gives none for the NUV, so Mg II has irispy's 53 mA.
  documented = [0.0318d, 0.0318d, 0.0286d, 0.0286d, 0.0286d, 0.053d]
  typical = 2

  openw, lun, 'iris_nonthermalwidth_formula.ecsv', /get_lun
  printf, lun, format='(A)', [ $
    '# %ECSV 1.0', '# ---', '# datatype:', '# - {name: ion, datatype: string}', $
    '# - {name: wavelength, unit: Angstrom, datatype: float64}', $
    '# - {name: instrumental_fwhm, unit: Angstrom, datatype: float64}', '# - {name: log_t, datatype: float64}', $
    '# - {name: fwhm, unit: Angstrom, datatype: float64}', '# - {name: velocity, unit: km / s, datatype: float64}', $
    '# meta: !!omap', $
    '# - comments:', $
    "#   - 'iris_nonthermalwidth(element, ion, wavelength, fwhm / (2 sqrt(ln 2)) / wavelength * c, instrumental_fwhm,'", $
    "#   - 'ti_use=log_t) (IDL 9.2, SSW 2026-09-28), with its c = 2.9979e5 km/s; velocity is -1 where fwhm is too narrow.'", $
    '# schema: astropy-2.0', 'ion wavelength instrumental_fwhm log_t fwhm velocity']
  ; iris_nonthermalwidth upper-cases its element and ion arguments, so they are passed as subscripts.
  for i = 0, 5 do begin
    velocity = iris_nonthermalwidth(elements[i], ions[i], wavelengths[i], fwhm[*, i] / (2 * sqrt(alog(2d))) $
      / wavelengths[i] * cc, instrumental[i], ti_use=log_t[i])
    foreach width, fwhm[*, i], j do printf, lun, '"' + elements[i] + ' ' + ions[i] + '"', shortest(wavelengths[i]), $
      shortest(instrumental[i]), shortest(log_t[i]), shortest(width), shortest(velocity[j]), format='(6(A, :, " "))'
  endfor
  free_lun, lun

  openw, lun, 'iris_nonthermalwidth_defaults.ecsv', /get_lun
  printf, lun, format='(A)', [ $
    '# %ECSV 1.0', '# ---', '# datatype:', '# - {name: ion, datatype: string}', $
    '# - {name: wavelength, unit: Angstrom, datatype: float64}', $
    '# - {name: instrumental_fwhm, unit: Angstrom, datatype: float64}', $
    '# - {name: fwhm, unit: Angstrom, datatype: float64}', '# - {name: log_t_max, datatype: float64}', $
    '# - {name: thermal_speed, unit: km / s, datatype: float64}', $
    '# - {name: velocity, unit: km / s, datatype: float64}', '# meta: !!omap', $
    '# - comments:', $
    "#   - 'iris_nonthermalwidth(element, ion, wavelength, fwhm / (2 sqrt(ln 2)) / wavelength * c, instrumental_fwhm,'", $
    "#   - 'ti_max=log_t_max, wt_v=thermal_speed) (IDL 9.2, SSW 2026-09-28): its default temperature and thermal speed,'", $
    "#   - 'and the instrumental FWHM of its header (none for the NUV, so 0.053 A for Mg II).'", $
    '# schema: astropy-2.0', 'ion wavelength instrumental_fwhm fwhm log_t_max thermal_speed velocity']
  for i = 0, 5 do begin
    width = fwhm[typical, i]
    velocity = iris_nonthermalwidth(elements[i], ions[i], wavelengths[i], width / (2 * sqrt(alog(2d))) / wavelengths[i] $
      * cc, documented[i], ti_max=log_t_max, wt_v=thermal_speed)
    printf, lun, '"' + elements[i] + ' ' + ions[i] + '"', shortest(wavelengths[i]), shortest(documented[i]), $
      shortest(width), shortest(log_t_max), shortest(thermal_speed), shortest(velocity), format='(7(A, :, " "))'
  endfor
  free_lun, lun
end
