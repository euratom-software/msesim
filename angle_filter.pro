function angle_filter, flt, lambda, angle, Neff

; Function calculates the effect (broadening) on the filtershape due to the angle spread of the fibre
; bundle.
;
; - input:  flt      : [nlambdax1] filtershape array (i.e. transmission curve)
;           lambda   : the wavelength vector in Angstrom
;           angle    : the (max.) half-angle in radians
;           Neff     : effective refractive index of the filter
;
; -returns: fltstruct: structure with 6 fields
;                      fltstruct.flt     = the corrected filter shape ([nlambdax1]-array)
;                      fltstruct.fltrc   = recentred, corrected filter shape ([nlambdax1]-array)
;                      fltstruct.wlshift = the shift of the centre in angstrom (always negative: blue shift)
;                      fltstruct.cwl     = new centre wavelength
;                      fltstruct.fwhm    = new FWHM
;                      fltstruct.T       = new peak transmission
;
;  v1.0 mdebock, 20/08/2007
;  v1.1 mdebock, 02/11/2007, wlshift added
;  v1.2 mdebock, 06/11/2007, renormalisation bug fixed
;  v2.0 mdebock, 17/10/2008, now the complete wavelength vector is an input. This lets us get rid of
;                            the hard coded wavelength of the transmission peak.

nlambda = n_elements(lambda)   ; size of the wavelength vector
dlambda = lambda[1]-lambda[0]  ; wavelength step

; we need the peak wavelength (centre of mass wavelength) of the filter.
cumfltold = total(flt,/cumulative)
oldcntidx = value_locate(cumfltold, 0.5*max(cumfltold))
if oldcntidx lt 0 then oldcntidx=0
if oldcntidx gt nlambda-1 then oldcntidx=nlambda-1
cwl = lambda[oldcntidx]

; loop through the incident angles (let's say in 100 steps)
steps=100
; the area emitting the light differs (get's bigger) when the angle is increased
area = 1.
; the new filter array
fltnew = fltarr(nlambda)
for i=1,steps do begin
  ; the angle of this ring
  theta=float(i)/steps*angle
  ; every wavelength of the is shifted by:
  wlshift = cwl * (1. - sqrt(1.-sin(theta)^2/Neff^2) )
  ; the wavelength shift corresponds with a shift in wavelength bins of
  binshift = round(wlshift/dlambda)
  if binshift lt nlambda then begin
    fltnew[0:(nlambda-1-binshift)] += area * flt[binshift:nlambda-1]
  endif
  ; get the area of the next ring
  area=float(i+1)^2-float(i)^2
endfor

; normalise new filter (there's no more light than that emitted the over the total area)
total_area=float(steps)^2
fltnew = fltnew/total_area


; find the shift of the peak
; we need the peak wavelength (centre of mass wavelength) of the filter.
cumfltnew = total(fltnew,/cumulative)
newcntidx = value_locate(cumfltnew, 0.5*max(cumfltnew))
if newcntidx lt 0 then newcntidx=0
if newcntidx gt nlambda-1 then newcntidx=nlambda-1
cwl = lambda[newcntidx]

; shift in wavelength bins
binshift = oldcntidx - newcntidx
if binshift lt 0       then message, 'Filter shifted to the red ??? NOT physical!!'
if binshift ge nlambda then message, 'Filter completely shifted off the wavelength axis! Increase size of the wavelength vector!'
; wavelength shift and new CWL
wlshift  = -dlambda*binshift
cwl      = lambda[newcntidx]

; new FWHM
mxflt = max(fltnew,mxidx)
if mxidx ne 0 then begin
  tmp   = min(abs(fltnew[0:(mxidx-1)] - 0.5*mxflt), halfidx1)
endif else halfidx1=0
tmp       = min(abs(fltnew[mxidx:(nlambda-1)] - 0.5*mxflt), halfidx2)
halfidx2 += mxidx
fwhm = lambda[halfidx2]-lambda[halfidx1]

; recentre fltnew
fltrc = fltarr(nlambda)
fltrc[binshift:nlambda-1] = fltnew[0:(nlambda-1-binshift)]


; return to total new filter and the recentred one
fltstruct = {flt:fltnew, fltrc:fltrc, wlshift:wlshift, cwl: cwl, fwhm:fwhm, T: mxflt}
return, fltstruct

end