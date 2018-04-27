function replace_nan, in, replace, infinity=infinity, nan=nan, sign=sign
;+
;  out = replace_nan(in, replace[, /infinity, /nan, sign=sign])
;
;  function replace all NaNs in the input variable with a (finite) replace value.
;
; :Params:
;    in      : in, required, type=(float,double, complex, string), scalar or array
;             The input variable. This can be a floating-point, double-precision, or complex
;             scalar or array expression. Strings are first converted to floating-point.
;             This function is meaningless for byte, integer, or longword arguments.
;    replace : in, required, type=(float,double, complex, string), scalar
;             The replace value. It should be of the same type as the input variable
; :Keywords:
;    infinity: in, optional, type=byte, scalar
;             Set this keyword to replace only infinite values
;    nan     : in, optional, type=byte, scalar
;             Set this keyword to replave only NaN values
;    sign    : in, optional, type=integer, scalar
;             If the infinity or nan keyword is set, then set this keyword will cause the function
;             only to replace positive infinities/NaNs is sign>0, only negative infinities/NaNs if
;             sign<0.
; :Returns:
;    out     : out, required, type=(float,double, complex, string), scalar or array
;             Returns a copy of the input variable with the NaNs replaced by the replace value.


  ; Default keyword values
  default, infinity, 0
  default, nan, 0
  default, sign, 0

  ; create the output variable
  out = in

  ; look for the non-finite values all non-finite values
  if ((infinity eq 0) && (nan eq 0)) || ((infinity eq 1) && (nan eq 1)) $
                                   then idx = where(~finite(in), cnt)
  ; look for infinity values only
  if (infinity eq 1) && (nan eq 0) then idx = where(finite(in, /infinity, sign=sign), cnt)
  ; look for NaN values only
  if (infinity eq 0) && (nan eq 1) then idx = where(finite(in, /nan, sign=sign), cnt)

  ; and replace them
  if cnt ne 0 then out[idx] = replace[0]

  return, out
end