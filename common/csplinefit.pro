; Non-linear LM-fit to a cubic spline function of an arbitrary number of nodes.
;------------------------------------------------------------------------------------------------------------
function cspline_to_fit, t_in, M_in
;+
;  out = cspline_to_fit(t_in, M_in)
;
;  Calculates the partial derivatives wrt the values at the M_in knot points (y_in) of a cubic spline at
;  the independend variable t_in. As the values at the knot points y_in are linearly independend, y_in
;  is not needed as an input parameter (we just need to know how many there are: M_in).
;  A common block is used for the x-coordinates of the spline knot points and - if requested - 
;  the gradients at the end points (in that case those end points are clamped rather than free).
;  The function is called by SVDFIT in order to fit this cubic spline to data points.
;
; :Params:
;     t_in : in, required, type=float, scalar or 1D array
;           independent variables
;     M_in : in, required, type=float, 1D array (at least size 2)
;           number of linearly independend fit parameters: i.e. the values y_in at the knot points
; :Returns:
;     out  : out, type=float, 2D array
;           the values of the partial derivatives to y_in at t_in (needed for fitting by SVDFIT)
;-

  ; x-coordinates of the nodes in common block
  common splinefit_common, xn_in, fit_idx, dx0_in, dxn_in

  ; reform the input parameters to x and y's
  x = xn_in
  y = fltarr(n_elements(xn_in))+1. ; a dummy set of knot point values
  t = t_in

  ; call the spline function with the 'fit'-keyword => returns also the partial derivatives to y.
  ; clamp the end points depending on whether dx0_in and dxn_in are finite numbers.
  if  finite(dx0_in) &&  finite(dxn_in) then out = cspline(x,y,t, dx0=dx0_in, dxn=dxn_in, /fit)
  if  finite(dx0_in) && ~finite(dxn_in) then out = cspline(x,y,t, dx0=dx0_in, /fit)
  if ~finite(dx0_in) &&  finite(dxn_in) then out = cspline(x,y,t, dxn=dxn_in, /fit)
  if ~finite(dx0_in) && ~finite(dxn_in) then out = cspline(x,y,t, /fit)

  ; we only need to return the partial derivatives with respect to y_in
  out = out[*,1:*]
  ; and from those only the once that to an y_in that is actually fitted
  out = out[*,fit_idx]

  ; return the result
  return, out

end

;------------------------------------------------------------------------------------------------------------

function csplinefit, xm, ym, xn, yn, errors=errors, dx0=dx0, dxn=dxn, yn_fix= yn_fix, double=double
;+
; out = csplinefit(xm, ym, xn, yn[, errors=errors, dx0=dx0, dxn=dxn, yn_fix= yn_fix, double=double])
;
; fits a spline with nodes [xn, yn] to the measurements [xm,ym]. yn is the initial guess at input and will
; be the fit result at output. The function returns the spline interpolation at xm.
;
; :Params:
;     xm    : in, required, type=float, 1D array
;            x-coordinates measurement
;     ym    : in, required, type=float, 1D array
;            y-coordinates measurement
;     xn    : in, required, type=float, 1D array
;            x-coordinates spline nodes
;     yn    : in/out, required, type=float, 1D array
;            y-coordinates spline nodes. Initial guess at input, fit result at output
; :Keywords:
;     errors: in, optional, type=float, 1D array
;            errors on ym
;     dx0   : in, optional, type=float, scalar
;            the derivative of S(x) at the first knot point: S'(x_in[0]) = dx0 <=> clamped spline knot.
;            If not set S''(x_in[0]) = 0 is used as boundary condition         <=> natural spline knot.
;     dxn   : in, optional, type=float, scalar
;            the derivative of S(x) at the last knot point : S'(x_in[n]) = dxn <=> clamped spline knot.
;            If not set S''(x_in[n]) = 0 is used as boundary condition         <=> natural spline knot.
;     yn_fix: in, optional, type=byte, 1D array
;            Array with the same number of elements as yn that sets the yn values that need to be kept fixed.
;     double: in, optional, type=byte, scalar
;            if set the calculation is done in double (higher precission)
; :Returns:
;     out  : the values of the fitted spline at xm
;-

  ; x-coordinates of the nodes in common block
  common splinefit_common, xn_in, fit_idx, dx0_in, dxn_in

  ; check the input data
  if n_elements(xm) ne n_elements(ym) then begin
    print, ''
    print, format='("ERROR: The number of x-coordinates (",I0,") is not equal to the number of number of y-coordinates (",'$
                 +'I0,") of the data to be fitted!")',n_elements(xm), n_elements(ym)
    print, ''
    return, -1.
  endif
  if keyword_set(errors) && (n_elements(xm) ne n_elements(errors)) then begin
    print, ''
    print, format='("ERROR: The number of errors on the data to be fitted (",I0,") is not equal to the number of data points (",'$
                 +'I0,")!")',n_elements(errors), n_elements(xm)
    print, ''
    return, -1.
  endif
  if n_elements(xn) ne n_elements(yn) then begin
    print, ''
    print, format='("ERROR: The number of x-coordinates (",I0,") is not equal to the number of number of y-coordinates (",'$
                 +'I0,") of the spline knot points!")',n_elements(xn), n_elements(yn)
    print, ''
    return, -1.
  endif

  ; get rid off NaN's in xm and in the errors and get rid of any negative or zero errors
  if keyword_set(errors) then begin
    idx = where( finite(xm) and finite(ym) and finite(errors) and (errors gt 0.), cnt)
    if cnt eq 0 then begin
      print, ''
      print, 'ERROR: No finite data points/errors were found or all errors are negative/zero!'
      print, ''
      return, -1.
    endif
    xm_in  = xm[idx]
    ym_in  = ym[idx]
    dym_in = errors[idx]
  endif else begin
    idx = where( finite(xm) and finite(ym), cnt)
    if cnt eq 0 then begin
      print, ''
      print, 'ERROR: No finite data points were found!'
      print, ''
      return, -1.
    endif
    xm_in  = xm[idx]
    ym_in  = ym[idx]
  endelse

  xn_in = xn
  yn_in = yn

  nn   = n_elements(xn_in)
  nm   = n_elements(xm_in)
  if nn ge nm then begin
    print, ''
    print, format='("ERROR: you are trying to fit ",I0," data points with a ",I0,"-node spline! ",A)',$
           nm, nn, 'The number of spline nodes should be at least 2 smaller than the number of data points!'
    print, ''
    return, -1.
  endif

  if n_elements(dx0) ne 0 then dx0_in = dx0 else dx0_in=!values.f_nan
  if n_elements(dxn) ne 0 then dxn_in = dxn else dxn_in=!values.f_nan

  ; y-coordinates of the nodes are the fit parameters
  if n_elements(yn_fix) ne n_elements(yn) then yn_fix=bytarr(n_elements(yn))
  fit_idx = where(~yn_fix, cnt)
  if cnt eq 0 then begin
    print, ''
    print, 'ERROR: all yn where requested to be fixed, hence no fitting can be done!'
    print, ''
    return, -1
  endif
  yn_fit = yn[fit_idx]

  ; Because the y-knots are linearly independend, we can use the fast and stable SVDfit algorithm
  ; However what SVDfit fits is: S_fit(x) = y_in[0] * dS(x)/dy_in[0] + ... + y_in[n] * dS[x]/dy_in[n]
  ; whereas the real spline is : S(x) = S_0(x) + y_in[0] * dS(x)/dy_in[0] + ... + y_in[n] * dS[x]/dy_in[n] = S_0(x) + S_fit(x)
  ; with S0(x) independent of y_in!
  ; So what we need to do is to find out S_0(xm), and subtract it from ym
  if  finite(dx0_in) &&  finite(dxn_in) then S = cspline(xn_in, yn_in, xm_in, dx0=dx0_in, dxn=dxn_in, /fit)
  if  finite(dx0_in) && ~finite(dxn_in) then S = cspline(xn_in, yn_in, xm_in, dx0=dx0_in, /fit)
  if ~finite(dx0_in) &&  finite(dxn_in) then S = cspline(xn_in, yn_in, xm_in, dxn=dxn_in, /fit)
  if ~finite(dx0_in) && ~finite(dxn_in) then S = cspline(xn_in, yn_in, xm_in, /fit)
  S_fit = fltarr(nm)
  for i=0,n_elements(fit_idx)-1 do S_fit += yn_in[fit_idx[i]] * S[*,(fit_idx[i]+1)]
  S_0   = S[*,0] - S_fit
  ym_in = ym_in - S_0

  ; do the fit
  if keyword_set(errors) then begin
    out = SVDfit( xm_in, ym_in, A=yn_fit, function_name='cspline_to_fit', measure_errors=dym_in, double=double, yfit=yfit)
  endif else begin
    out = SVDfit( xm_in, ym_in, A=yn_fit, function_name='cspline_to_fit', double=double, yfit=yfit)
  endelse

  ; return the result
  yn[fit_idx] = out
  yfit = yfit + S_0
  return, yfit
end