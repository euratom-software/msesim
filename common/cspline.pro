function cspline, x_in, y_in, x_out, dx0=dx0, dxn=dxn, fit=fit, coeff=coeff
;+
; out = cspline(x_in, y_in, x_out [, da=da, db=db, fit=fit])
;
; Free and clampled cubic spline function based on the algorithm described in cspline.pdf [~/my_idl/common/cspline.pdf]
; It performs a free (natural) or  clamped,cubic spline interpolation form the (x_in,y_in)-nodes to the new x-axis x_out.
; If the 'fit'-keyword is set, the function will also return the partial derivatives to y_in which are
; needed for the fitting. If the 'fit'-keyword is not set than only the interpolated spline data is returned.
; The 'dx0' and 'dxn'-keywords give the values of the derivatives at the end points, in case of a clamped
; spline: S'(x_in[0])=dx0 and/or S'(x_in[n])=dxn. If omitted the free or natural boundary condition is used:
; S''(x_in[0])=S''(x_in[n])=0. It is possible to only clamp the spline at one end point.
;
;   S(x) = Si(x) = S(x) = S[i](x) = a[i]+ b[i]*(x - x_in[i]) + c[i]*(x - x_in[i])^2 + d[i]*(x - x_in[i])^3
;     for x_in[i] <= x <= x_in[i+1] and i=0,...n-1.
;     with
;      (1) S[i](x_in[i])       = y[i]
;      (2) S[i+1](x_in[i+1])   = S[i](x_in[i+1])
;      (3) S'[i+1](x_in[i+1])  = S'[i](x_in[i+1])
;      (4) S''[i+1](x_in[i+1]) = S''[i](x_in[i+1])
;      (5) S''(x_in[0]) = S''(x_in[n]) = 0           => natural spline <=> boundary keyword not set
;         or
;          S'(x_in[0]) = dx0 and S'(x_in[n]) = dxn   => clamped spline <=> boundary keyword set to end derivatives
;
; :Params:
;     x_in : in, required, type=float, 1D array (at least size 2)
;           x-coordinates of the spline knots
;     y_in : in, required, type=float, 1D array (at least size 2)
;           y-coordinates of the spline knots
;     x_out: in, required, type=float, scalar or 1D array
;           new x-coordinates
; :Keywords:
;     dx0  : in, optional, type=float, scalar
;           the derivative of S(x) at the first knot point: S'(x_in[0]) = dx0 <=> clamped spline knot.
;           If not set S''(x_in[0]) = 0 is used as boundary condition         <=> natural spline knot.
;     dxn  : in, optional, type=float, scalar
;           the derivative of S(x) at the last knot point : S'(x_in[n]) = dxn <=> clamped spline knot.
;           If not set S''(x_in[n]) = 0 is used as boundary condition         <=> natural spline knot.
;     fit  : in, optional, type=byte, scalar
;           if set the function will also return the partial derivatives to y
; :Returns:
;     out  : out, type=float, 1D array (if keyword fit not set) or 2D array (if keyword fit set)
;           if keyword fit not set - the values of the spline interpolation. 1D array with the same size as x_out
;           if keyword fit set     - 2D array [n_elements(x_out), n_elements(y_in)+1].
;                                    first row contains the values of the spline interpolation. The following
;                                    n_elements(y_in)-rows contain the partial derivatives to each y_in.
;-

  ; number of cubic polynomials in the spline
  n = n_elements(x_in)-1

  ; The a-coefficients are easily found from condition (1): S[i](x_in[i])=y_in[i] <=> a[i]=y_in[i]
  a = y_in

  ; initialise the b-, c- and d-coefficients
  b = fltarr(n+1)
  c = fltarr(n+1)
  d = fltarr(n+1)

  ; Using the conditions (2),(3) and (4) you find following expression for c[i] (see equation (24) in cspline.pdf):
  ;  r[i] * c[i-1] + s[i]*c[i] + t[i]*c[i+1] = m[i], for i=1,...,n-1
  ; with:
  ;  h[i] = x_in[i+1]-x_in[i]
  ;  r[i] = h[i-1]
  ;  s[i] = 2.*(h[i]+h[i-1])
  ;  t[i] = h[i]
  ;  m[i] = 3./h[i]*(a[i+1]-a[i]) - 3./h[i-1]*(a[i]-a[i-1])
  h        = x_in[1:n]-x_in[0:n-1]
  m        = fltarr(n+1)
  r        = fltarr(n+1)
  s        = fltarr(n+1)
  t        = fltarr(n+1)
  if n gt 1 then begin
    r[1:n-1] = h[0:n-2]
    s[1:n-1] = 2.*(h[1:n-1]+h[0:n-2])
    t[1:n-1] = h[1:n-1]
    m[1:n-1] = 3./h[1:n-1]*(a[2:n]-a[1:n-1]) - 3./h[0:n-2]*(a[1:n-1]-a[0:n-2])
  endif

  ; from this b[i] and d[i] can be calculated using equations (18) and (20) in cspline.pdf:
  ;  d[i] = (c[i+1]-c[i])/(3.*h[i])                       , i=0,...,n-1
  ;  b[i] = (a[i+1]-a[i])/h[i] - h[i]/3.*(2.*c[i]+c[i+1]) , i=0,...,n-1

  ; to find c[0] and c[n] we need to use the boundary conditions, which will return
  ;   s[0], t[0], m[0] for the equation that returns c[0] : s[0]*c[0]   + t[0]*c[1] = m[0]
  ;   r[n], s[n], m[n] for the equation that returns c[n] : r[n]*c[n-1] + s[n]*c[n] = m[n]
  if n_elements(dx0) eq 0 then begin
    s[0] = 1.
    t[0] = 0.
    m[0] = 0.
  endif else begin
    s[0] = 2.*h[0]
    t[0] = h[0]
    m[0] = 3./h[0] * (a[1]-a[0]) - 3.*dx0
  endelse

  if n_elements(dxn) eq 0 then begin
    r[n] = 0.
    s[n] = 1.
    m[n] = 0.
  endif else begin
    r[n] = h[n-1]
    s[n] = 2.*h[n-1]
    m[n] = 3.*dxn - 3./h[n-1] * (a[n]-a[n-1])
  endelse

  ; We can now solve the set of equations: r[i]*c[i-1] + s[i]*c[i] + t[i]*c[i+1] = m[i], for i=0,...,n
  ; with a tridiagonal matrix algorithm. Recursively substituting
  ;  u[i] = u[i-1]*s[i]-r[i]*v[i-1] and u[0]= s[0]
  ;  v[i] = u[i-1]*t[i]             and v[0]= t[0]
  ;  w[i] = u[i-1]*m[i]-r[i]*w[i-1] and w[0]= m[0]
  ; the set of equations reduces to      : u[i]*c[i]+v[i]*c[i+1] =w[i]                 , for i=0,...,n
  u = fltarr(n+1)
  v = fltarr(n+1)
  w = fltarr(n+1)
  u[0] = s[0]
  v[0] = t[0]
  w[0] = m[0]
  for i=1,n do begin
    u[i] = u[i-1]*s[i]-r[i]*v[i-1]
    v[i] = u[i-1]*t[i]
    w[i] = u[i-1]*m[i]-r[i]*w[i-1]
  endfor

  ; because t[n] is by definition 0, v[n] =0 as well and the solution for c[n] is trivial:
  c[n] = w[n]/u[n]
  ; backwards substitution then yield every c[i]: c[i]=(w[i]-v[i]*c[i+1])/u[i], and also b[i] and d[i]
  for i=n-1,0,-1 do begin
    c[i] = (w[i]-v[i]*c[i+1])/u[i]
    d[i] = (c[i+1]-c[i])/(3.*h[i])
    b[i] = (a[i+1]-a[i])/h[i] - h[i]/3.*(2.*c[i]+c[i+1])
  endfor

  ; we now have all the coefficients so we can calculate the spline
  ; calculate spline
  y_out = fltarr(n_elements(x_out))
  for i=0,n-1 do begin
    ; find the interval for x_out
    idx = where((x_out ge x_in[i]) and (x_out le x_in[i+1]),cnt)
    if i eq 0   then idx = where(x_out le x_in[i+1],cnt)
    if i eq n-1 then idx = where(x_out ge x_in[i]  ,cnt)
    if cnt ne 0 then begin
      y_out[idx] = a[i]+b[i]*(x_out[idx]-x_in[i])+c[i]*(x_out[idx]-x_in[i])^2+d[i]*(x_out[idx]-x_in[i])^3
    endif
  endfor

  ; stop here if the 'fit'-keyword is not set
  if ~keyword_set(fit) then return, y_out


  ; if the 'fit'-keyword is set we need to calculate the polynomials of the partial derivatives to y_in
  da = fltarr(n,n+1)
  db = fltarr(n,n+1)
  dc = fltarr(n,n+1)
  dd = fltarr(n,n+1)
  for i=0,n do begin
    ; da/dy is simple, because a=y
    da_dy    = fltarr(n+1)
    da_dy[i] = 1.

    ; secondly get dm/dy
    dm_dy    = fltarr(n+1)
    for j=1,n-1 do begin
      dm_dy[j] = 3./h[j] * (da_dy[j+1]-da_dy[j]) - 3./h[j-1] * (da_dy[j]-da_dy[j-1])
    endfor
    if n_elements(dx0) eq 0 then dm_dy[0] = 0. else begin
      dm_dy[0] =  3./h[0] * (da_dy[1]-da_dy[0])
    endelse
    if n_elements(dxn) eq 0 then dm_dy[n] = 0. else begin
      dm_dy[n] = -3./h[n-1] * (da_dy[n]-da_dy[n-1])
    endelse

    ; then get dw/dy
    dw_dy    = fltarr(n+1)
    dw_dy[0] = dm_dy[0]
    for j=1,n do begin
      dw_dy[j] = u[j-1]*dm_dy[j]-r[j]*dw_dy[j-1]
    endfor

    ; now we can get dc/dy, db/dy and dd/dy
    dc_dy    = fltarr(n+1)
    db_dy    = fltarr(n+1)
    dd_dy    = fltarr(n+1)
    dc_dy[n] = dw_dy[n]/u[n]
    for j=n-1,0,-1 do begin
      dc_dy[j] = (dw_dy[j]-v[j]*dc_dy[j+1])/u[j]
      dd_dy[j] = (dc_dy[j+1]-dc_dy[j])/(3.*h[j])
      db_dy[j] = (da_dy[j+1]-da_dy[j])/h[j] - h[j]/3.*(2.*dc_dy[j]+dc_dy[j+1])
    endfor

    ; finally store it all
    da[*,i] = da_dy[0:n-1]
    db[*,i] = db_dy[0:n-1]
    dc[*,i] = dc_dy[0:n-1]
    dd[*,i] = dd_dy[0:n-1]
  endfor

  ; calculate splines of the partial derivatives to y_in
  dy_out = fltarr(n_elements(x_out),n+1)
  for i=0,n-1 do begin
    ; find the interval for x_out
    idx = where((x_out ge x_in[i]) and (x_out le x_in[i+1]),cnt)
    if i eq 0   then idx = where(x_out le x_in[i+1],cnt)
    if i eq n-1 then idx = where(x_out ge x_in[i]  ,cnt)
    if cnt ne 0 then begin
      for j=0,n do begin
        dy_out[idx,j] = da[i,j]+db[i,j]*(x_out[idx]-x_in[i])+dc[i,j]*(x_out[idx]-x_in[i])^2+dd[i,j]*(x_out[idx]-x_in[i])^3
      endfor
    endif
  endfor



  ; return the output
  out            = fltarr(n_elements(x_out),n+2)
  out[*,0]       = y_out
  out[*,1:(n+1)] = dy_out
  return, out

end