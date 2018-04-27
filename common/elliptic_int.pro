function elliptic_int, p, k, incomplete=incomplete
;+
;  out = elliptic_int(p, k, incomplete=incomplete)
;
;     Function: ELLIPTIC_INT
;     Version : 1.00
;     Author  : M. De Bock
;     Date    : October - 2009
;
; A function to calculate the first and second order elliptic integral.
; It uses the standard tabulated integration function of IDL: INT_TABULATED.
;
; If the 'incomplete' keyword is omitted or set to zero the complete integral is calculated.
; If the 'incomplete' keyword is set to t, then the incomplete elliptical integral for t is calculated.
; The complete intergral is that at t=pi/2:
;     elliptic_int(p, k) == elliptic_int(p, k, incomplete=!pi/2.)
;
; The formula's for the elliptical integrals of 1st and second kind are:
;
;                                                  t
;                                                 /
;   1st kind:  elliptic_int(p, 1, incomplete=t) = |  1/sqrt( 1 - p^2 * sin(theta)^2 ) dtheta
;                                                 /
;                                                 0
;
;                                                  t
;                                                 /
;   2nd kind:  elliptic_int(p, 2, incomplete=t) = |  sqrt( 1 - p^2 * sin(theta)^2 ) dtheta
;                                                 /
;                                                 0
;
;
; The complete integrals were checked againts R. Martin's ELLIPTIC routine
;
; :Params:
;    p          : in, required, type=float/double
;                The elliptic modulus. It has to satisfy 0 < p^2 < 1
;    k          : in, required, type=byte
;                The kind of elliptic integral: should be either 1 or 2.
; :Keywords:
;    incomplete : in, type=float/double
;                if not set or set to zero, then the complete elliptical integral is calculated.
;                else the incomplete elliptical integral for t=incomplete is calculated.
;                rather than the WinSpec wavelength calibration (is more accurate!)
; :Returns:
;    out        : the function value
;-

  ; Check input parameters
  show_usage=0
  if n_params() ne 2 then show_usage=1 $
  else begin
    if (p lt 0.) || (p gt 1.)   then show_usage=1
    if ~((k eq 1) || (k eq 2)) then show_usage=1
  endelse
  if show_usage then begin
    print, ''
    print, '  Function: ELLIPTIC_INT'
    print, '  Version : 1.00'
    print, '  Author  : M. De Bock'
    print, '  Date    : October - 2009'
    print, ''
    print, '  A function to calculate the first and second order elliptic integral.'
    print, ''
    print, '  Usage:'
    print, '    out = elliptic_int(p, k, incomplete=t)'
    print, ''
    print, '       - p            : elliptic modulus (0 < p < 1) [required]'
    print, '       - k            : order (1 or 2)               [required]'
    print, '       - incomplete=t : if not set or set to zero, then the complete elliptical integral is returned.'
    print, '                         else the incomplete elliptical integral for t=incomplete is returned.'
    print, ''
    return, !values.f_NaN
  endif

  ; check incomplete keyword and find out whether we need to return a float rather than a double
  dofloat = 0
  pf = (size(p,/type) eq 4)
  if keyword_set(incomplete) then begin
    t  = incomplete
    tf = (size(t,/type) eq 4)
    if pf && tf then dofloat = 1
  endif else begin
    if pf then begin
      dofloat = 1
      t = !pi/2.
    endif else t = 3.14159265358979323846D / 2.0D
  endelse

  ; define the theta array. To get an accurate calculation we do everything in double, and we
  ; use an extensive theta array for performing the integration
  dtheta = !pi/1e4
  ntheta = floor(t/dtheta)
  theta  = dindgen(ntheta)/(ntheta-1.) * t

  ; Elliptical integral of the first kind
  if k eq 1 then begin
    int = int_tabulated(theta, 1./sqrt(1.-p^2 * sin(theta)^2 ), /double)
    if dofloat then int=float(int)                          ; convert back from double to float is p is a float
  endif

  ; Elliptical integral of the second kind
  if k eq 2 then begin
    int = int_tabulated(theta, sqrt(1.-p^2 * sin(theta)^2 ), /double)
    if dofloat then int=float(int)                          ; convert back from double to float is p is a float
  endif

  return, int
end