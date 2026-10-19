; Writes the fit_iris_xput outputs read by test_fit_xput_lite_idl in irispy/utils/tests/test_response.py.
; Needs SolarSoft with this path:
;   !path = expand_path('+~/ssw/iris/idl') + ':' + expand_path('+~/ssw/gen/idl') + ':' + !path
; Run it in a directory holding fit_iris_xput_input.sav: tt, c_f_time (cftime) and coeffs_fuv of
; response version 9 (iris_sra_c_20231106.geny), saved by iris_get_response for 2025-08-05.
;
; Each file holds rr = fltarr(1, 3), the FUV throughput at the three wavelengths of c_f_lambda,
; computed as iris_get_response does (IDL 9.2, SolarSoft of 2026-09-28).
; tt (2025-08-05 20:44:06 UTC) is in the last calibration interval and 2013-10-22 21:00 in the
; first. 2020-11-15 00:00 is in the gap between the intervals 2018-06-15 to 2020-10-01 and
; 2021-01-01 to 2021-11-28, where the fit blends the two.

pro make_fit_iris_xput_output
  compile_opt idl2
  restore, 'fit_iris_xput_input.sav'
  times = list(tt, '2013-10-22T21:00', '2020-11-15T00:00')
  files = ['fit_iris_xput_output.sav', 'fit_iris_xput_output_20131022.sav', 'fit_iris_xput_output_20201115.sav']
  foreach time, times, i do begin
    rr = fltarr(1, 3)
    for j = 0, 2 do rr[*, j] = fit_iris_xput(time, cftime, coeffs_fuv[*, *, j])
    save, rr, filename=files[i]
  endforeach
end
