; Non-linear LM-fit to a cubic spline function of an arbitrary number of nodes.


function spline_mdb, x_in, y_in, t_in, fit=fit
;+
; out = spline_mdb(x_in, y_in, t_in)
;
; my own cubic spline function based on the algorithm described in wikipedia:
; [http://en.wikipedia.org/wiki/Spline_(mathematics)#General_Expression_For_a_C2_Interpolating_Cubic_Spline]
; It performs a clamped, cubic spline interpolation form the (x_in,y_in)-nodes to the new x-axis t_in.
; If the 'fit'-keyword is set, the function will also return the partial derivatives to y which are 
; needed for the fitting. If the 'fit'-keyword is not set than only the interpolate spline data is returned.
;
; :Params:
;     x_in : in, required, type=float, 1D array (at least size 2)
;           x-coordinates of the spline nodes
;     y_in : in, required, type=float, 1D array (at least size 2)
;           y-coordinates of the spline nodes
;     t_in : in, required, type=float, scalar or 1D array
;           new x-coordinates
; :Keywords:
;     fit  : in, type=byte, scalar
;           if set the function will also return the partial derivatives to y
; :Returns:
;     out  : out, type=float, 1D array (if keyword fit not set) or 2D array (if keyword fit set)
;           if keyword fit not set - the values of the spline interpolation. 1D array with the same size as t_in
;           if keyword fit set     - 2D array [n_elements(t_in), n_elements(y_in)+1].
;                                    first row contains the values of the spline interpolation. The following
;                                    n_elements(y_in)-rows contain the partial derivatives to each y_in.
;-

  ; reform the input parameters to x and y's
  x = x_in
  y = y_in
  t = t_in
  n = n_elements(x)-1

  ; calculate spline polynomials
  a = y
  b = fltarr(n)
  d = fltarr(n)
  h = x[1:n]-x[0:n-1]

  alpha = fltarr(n+1)
  alpha[1:n-1] = 3./h[1:n-1] * (a[2:n]-a[1:n-1]) - 3./h[0:n-2] * (a[1:n-1]-a[0:n-2])

  c  = fltarr(n+1)
  l  = fltarr(n+1)
  mu = fltarr(n+1)
  z  = fltarr(n+1)

  l[0]  = 1.
  mu[0] = 0.
  z[0]  = 0.
  for i=1,n-1 do begin
    l[i]  = 2.*(x[i+1]-x[i-1])-h[i-1]*mu[i-1]
    mu[i] = h[i]/l[i]
    z[i]  = (alpha[i]-h[i-1]*z[i-1])/l[i]
  endfor
  l[n]  = 1.
  c[n]  = 0.
  z[n]  = 0.
  for j=n-1,0,-1 do begin
    c[j] = z[j]-mu[j]*c[j+1]
    b[j] = (a[j+1]-a[j])/h[j] - h[j]*(c[j+1]+2.*c[j])/3.
    d[j] = (c[j+1]-c[j])/(3.*h[j])
  endfor
  S = fltarr(n,4)
  S[*,0] = a[0:n-1]
  S[*,1] = b[0:n-1]
  S[*,2] = c[0:n-1]
  S[*,3] = d[0:n-1]

  ; calculate spline
  y_out = fltarr(n_elements(t))
  for i=0,n-1 do begin
    ; find the interval for t
    idx = where((t ge x[i]) and (t lt x[i+1]),cnt)
    if i eq 0 then idx = where(t lt x[i+1],cnt)
    if i eq n-1 then idx = where(t ge x[i],cnt)
    if cnt ne 0 then begin
      y_out[idx,0] = S[i,0]+S[i,1]*(t[idx]-x[i])+S[i,2]*(t[idx]-x[i])^2+S[i,3]*(t[idx]-x[i])^3
    endif
  endfor

  ; stop here if the 'fit'-keyword is not set
  if ~keyword_set(fit) then return, y_out


  ; if the 'fit'-keyword is set we need to calculate the polynomials of the partial derivatives to y
  dS = fltarr(n,4,n+1)
  for i=0,n do begin
    ; da/dy is simple, because a=y
    da_dy    = fltarr(n+1)
    da_dy[i] = 1.
    ; secondly get dalpha/dy
    dalpha_dy = fltarr(n+1)
    for j=1,n-1 do begin
      dalpha_dy[j] = 3./h[j] * (da_dy[j+1]-da_dy[j]) - 3./h[j-1] * (da_dy[j]-da_dy[j-1])
    endfor
    ; then get dz/dy
    dz_dy    = fltarr(n+1)
    dz_dy[0] = 0.
    dz_dy[n] = 0.
    for j=1,n-1 do begin
      dz_dy[j] = (dalpha_dy[j]-h[j-1]*dz_dy[j-1])/l[j]
    endfor
    ; now we can get dc/dy, db/dy and dd/dy
    dc_dy    = fltarr(n+1)
    db_dy    = fltarr(n+1)
    dd_dy    = fltarr(n+1)
    dc_dy[n] = 0
    for j=n-1,0,-1 do begin
      dc_dy[j] = dz_dy[j]-mu[j]*dc_dy[j+1]
      db_dy[j] = (da_dy[j+1]-da_dy[j])/h[j] - h[j]*(dc_dy[j+1]+2.*dc_dy[j])/3.
      dd_dy[j] = (dc_dy[j+1]-dc_dy[j])/(3.*h[j])
    endfor
    ; finally store it all
    dS[*,0,i] = da_dy[0:n-1]
    dS[*,1,i] = db_dy[0:n-1]
    dS[*,2,i] = dc_dy[0:n-1]
    dS[*,3,i] = dd_dy[0:n-1]
  endfor


  ; calculate splines
  dy_out = fltarr(n_elements(t),n+1)
  for i=0,n-1 do begin
    ; find the interval for t
    idx = where((t ge x[i]) and (t lt x[i+1]),cnt)
    if i eq 0 then idx = where(t lt x[i+1],cnt)
    if i eq n-1 then idx = where(t ge x[i],cnt)
    if cnt ne 0 then begin
      for j=0,n do begin
        dy_out[idx,j] = dS[i,0,j]+dS[i,1,j]*(t[idx]-x[i])+dS[i,2,j]*(t[idx]-x[i])^2+dS[i,3,j]*(t[idx]-x[i])^3
      endfor
    endif
  endfor
  ; return the output
  out            = fltarr(n_elements(t),n+2)
  out[*,0]       = y_out
  out[*,1:(n+1)] = dy_out
  return, out

end

;------------------------------------------------------------------------------------------------------------

function spline_to_fit, t_in, A_in
;+
;  out = spline_to_fit(t_in, A_in)
;
;  Calculates a clamped cubic spline based on the independend variable t_in and the parameters A_in.
;  These A_in parameters are th y-coordinates of the spline-nodes. A common block is used for the
;  x-coordinates of the spline-nodes. The function is called by LMFIT in order to fit this cubic spline
;  to data points.
;
; :Params:
;     t_in : in, required, type=float, scalar or 1D array
;           independent variables
;     A_in : in, required, type=float, 1D array (at least size 2)
;           fit parameters
; :Returns:
;     out  : out, type=float, 2D array
;           the value of the spline at t_in, and the partial derivatives to y needed for the fitting
;-

  ; x-coordinates of the nodes in common block
  common splinefit_common, xn_splinefit

  ; reform the input parameters to x and y's
  x = xn_splinefit
  y = A_in
  t = t_in

  ; call the spline function with the 'fit'-keyword => returns also the partial derivatives to y
  out = spline_mdb(x,y,t, /fit)
  ; return the result
  return, out

end

;------------------------------------------------------------------------------------------------------------

function splinefit, xm, ym, xn, yn, errors=errors, double=double
;+
; out = splinefit(xm, ym, xn, yn, errors=errors)
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
;     double: in, optional, type=byte, scalar
;            if set the calculation is done in double (higher precission)
; :Returns:
;     out  : the values of the fitted spline at xm
;-

  nn   = n_elements(xn)
  nm   = n_elements(xm)
  if nn ge nm then begin
    print, format='("ERROR: you are trying to fit ",I0," data points with a ",I0,"-node spline! ",A)',$
           nm, nn, 'The number of spline nodes should be at least 2 smaller than the number of data points!'
  endif

  ; x-coordinates of the nodes in common block
  common splinefit_common, xn_splinefit
  xn_splinefit = xn

  ; y-coordinates of the nodes are the fit parameters
  A    = yn

  ; do the fit
  out = LMfit( xm, ym, A, function_name='spline_to_fit', measure_errors=errors, double=double)

  ; return the result
  yn = A
  return, out
end