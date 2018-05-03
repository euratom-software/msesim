;; Copyright 2018 Euratom/CCFE/TuE

;; Permission is hereby granted, free of charge, to any person obtaining a copy of this
;; software and associated documentation files (the "Software"), to deal in the Software
;; without restriction, including without limitation the rights to use, copy, modify,
;; merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
;; permit persons to whom the Software is furnished to do so, subject to the following
;; conditions:

;; The above copyright notice and this permission notice shall be included in all copies
;; or substantial portions of the Software.

;; THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
;; INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
;; PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
;; HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
;; CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
;; OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

function calc_equi, pos, R0, a, shafr, elong, Bphi, q0, qa, qindex, Bp0, Bpa, Bpindex

; this function calculates the B-field ([3x1] vector Bfld), the flux surface index (psi) and the poloidal
; angle (theta) for a given point pos and the given tokamak and field parameters (R0, q0, etc.). It also
; returns the magnetic axis. 
;
; - input:  pos     [3xn] vector-array containing the xyz-coordinates of the point
;           R0      major radius of the geometrical centre of the tokamak (also centre LCFS)
;           a       minor radius (i.e. radius LCFS) 
;           shafr   Shafranov shift (i.e. magnetic axis at R0+shafr)
;           elong   elongation (i.e. height of the plasma is 2*a*elong)
;           Bphi    Toroidal field at R0
;           q0      q at the magnetic axis
;           qa      q at the plasma edge (LCFS)
;           qindex  q = q0 + (qa-q0)*sqrt(psi)^qindex
;           Bp0     (Toroidal) paramagnetic field at the magnetic axis
;           Bpa     (Toroidal) paramagnetic field at the plasma edge (LCFS)
;           Bpindex Bp = Bp0 + (Bpa-Bp0)*sqrt(psi)^Bpindex
;
; -returns: equi    a structure with fields: Bfld  a [3xn] vector-array with the xyz-components of the B-field at position 'pos'
;                                            psi   a [nx1] array with values in between 0 (magn. axis) and 1 (LCFS) that
;                                                  indicates on which flux surface 'pos' lies
;                                            theta a [nx1] array with values in between 0 and 2*pi that indicates the poloidal angle
;                                            Rm    a scalar that contains the major radius of the magnetic axis
;  v1.0 mdebock, 05/06/2007
;
;  calculate the equilibrium at all the gridpoints at once, instead of calculating the B-field for just one grid point
;  (so this calculation can be taken out of the loop!)
;  v1.1 mdebock, 22/06/2007


; first transform pos into major radius R and Z
R = sqrt(pos[0,*]^2 + pos[1,*]^2)
Z = pos[2,*]

; initialise the sqrt(psi) and theta vectors.
; REMARK: we call the variable 'psi', but at this point it actually is sqrt(psi), 
;         which is usuallt approx. the normalised minor radius and in this case exactly 
;         the normalised minor radius.
n     = n_elements(R)
psi   = fltarr(1,n) + 1.0
theta = fltarr(1,n)

; if psi and theta are known, than the corresponding R and Z are, 
; assuming the Shafranov shift goes quadratically to 0 when psi goes to 1:
;  R = R0 + shafr*(1 - psi^2) + psi*a*cos(theta)	(1)
;  Z = psi*a*elong*sin(theta)				(2)

; In the rare case that (R,Z) is exactly on the magnetic axis (R=R0+shafr,Z=0)
; psi and theta are 0.
; REMARK: we're working with floats that are never exactly 0, so we'll work with a tolerance EPS
EPS =1e-5	; i.e. distance less than 0.01 mm)
idx = where( (abs(R-R0-shafr) le EPS) AND (abs(Z) le EPS),  count)
if (count ne 0) then begin
  psi[idx]   = 0.0
  theta[idx] = 0.0
end


; if Z=0 and R>R0+shafr then theta=0
idx = where( ((R-R0-shafr) gt EPS) AND (abs(Z) le EPS),  count)
if (count ne 0) then begin
  theta[idx] = 0.0
  ; R = R0 + shafr*(1 - psi^2) + psi*a
  ; shafr*psi^2 - a*psi + [R - shafr - R0] = 0
  Rtmp       = transpose(R[idx])
  coeff      = [Rtmp - shafr - R0,      $
                replicate(-a,1,count),  $
                replicate(shafr,1,count)]
  psitmp     = real_part(quadroots(coeff))

  ; because psi lies in the range [0,1]
  psiidx     = where((psitmp le 1+EPS) AND (psitmp ge 0-EPS),cnt)
  if (cnt ne 0) then begin
    if (size(psitmp,/n_dimensions) eq 1) then begin	; if there is just one point in idx, then psitmp is [2x1]
      psi[idx] = psitmp[psiidx]				; and the function array_indices doesn't work properly
    endif else begin
      psiidx = ARRAY_INDICES(psitmp,psiidx)
      psi[idx[psiidx[1,*]]] = psitmp[psiidx[0,*],psiidx[1,*]]
    endelse
  endif
endif

; if Z=0 and R<R0+shafr then theta=pi
idx = where( ((R-R0-shafr) lt -EPS) AND (abs(Z) le EPS),  count)
if (count ne 0) then begin
  theta[idx] = !pi
  ; R = R0 + shafr*(1 - psi^2) - psi*a
  ; shafr*psi^2 + a*psi + [R - shafr - R0] = 0
  Rtmp       = transpose(R[idx])
  coeff      = [Rtmp - shafr - R0,      $
                replicate(a,1,count),   $
                replicate(shafr,1,count)]
  psitmp     = real_part(quadroots(coeff))

  ; because psi lies in the range [0,1]
  psiidx     = where((psitmp le 1+EPS) AND (psitmp ge 0-EPS),cnt)
  if (cnt ne 0) then begin
    if (size(psitmp,/n_dimensions) eq 1) then begin	; if there is just one point in idx, then psitmp is [2x1]
      psi[idx] = psitmp[psiidx]				; and the function array_indices doesn't work properly
    endif else begin
      psiidx = ARRAY_INDICES(psitmp,psiidx)
      psi[idx[psiidx[1,*]]] = psitmp[psiidx[0,*],psiidx[1,*]]
    endelse
  endif
endif

; In all the other cases (ie. Z <> 0) can solve the set of equations
; by substituting (2) into (1):
idx = where(abs(Z) gt EPS,  count)
if (count ne 0) then begin
  ;  R = R0 + shafr*[1 - Z^2/(a*elong)^2 * (1 + x^2) ] + Z/elong * x
  ;  [shafr*Z^2/(a*elong)^2]*x^2 - Z/elong * x + [R - shafr + shafr*Z^2/(a*elong)^2 - R0] = 0
  ; with x=cotan(theta)
  Rtmp       = transpose(R[idx])
  Ztmp       = transpose(Z[idx])
  coeff      = [Rtmp - shafr + shafr*Ztmp^2/(a*elong)^2 - R0,$
                -Ztmp/elong,                                   $
                shafr*Ztmp^2/(a*elong)^2                       ]
  x          = real_part(quadroots(coeff))

  thetatmp = fltarr(2,count)
  thetaidx = where(x eq 0.0,cnt)
  if (cnt ne 0) then begin
    thetatmp[thetaidx] = !pi/2
  endif
  thetaidx = where(x ne 0.0,cnt)
  if (cnt ne 0) then begin
    thetatmp[thetaidx] = atan(1/x[thetaidx])
  endif
  thetatmp = [thetatmp, thetatmp+!pi]	; there are 4 possible angles theta, because each atan(x) that 2 solutions
  thetaidx = where(thetatmp lt 0.0,cnt)
  if (cnt ne 0) then begin
    thetatmp[thetaidx] = thetatmp[thetaidx] + 2*!pi	; make sure theta lies in the range [0,2*pi] (and not [-pi/2,3*pi/2])
  endif

  ; the corresponding values of psi are: 
  psitmp   = rebin(Ztmp,4,count)/(a*elong*sin(thetatmp))

  ; because psi lies in the range [0,1]
  psiidx     = where((psitmp le 1+EPS) AND (psitmp ge 0-EPS),cnt)
  if (cnt ne 0) then begin
    if (size(psitmp,/n_dimensions) eq 1) then begin	; if there is just one point in idx, then psitmp is [2x1]
      psi[idx]   = psitmp[psiidx]			; and the function array_indices doesn't work properly
      theta[idx] = thetatmp[psiidx]
    endif else begin
      psiidx = ARRAY_INDICES(psitmp,psiidx)
      psi[idx[psiidx[1,*]]] = psitmp[psiidx[0,*],psiidx[1,*]]
      theta[idx[psiidx[1,*]]] = thetatmp[psiidx[0,*],psiidx[1,*]]
    endelse
  endif
endif



; from psi we find q
q = q0 + (qa-q0)*psi^qindex

; and from q we find the poloidal B-field Btheta
Rpsi   = R0+shafr*(1-psi^2)		; major radius of the centre of the flux surface
Btheta = a*psi*(Bphi*R0/Rpsi)/(Rpsi*q)

; the toroidal B-field Bphi at these positions are
BphiR  = Bphi*R0/R			; vacuum field
Bpara  = Bp0 + (Bpa-Bp0)*psi^Bpindex	; paramagnetic field
BphiR += Bpara				; total toroidal field

; from BphiR, Btheta, theta and phi we can derive the B-field vector
phi        = atan(pos[1,*],pos[0,*])
Bphi_vec   = rebin(BphiR,3,n)  * [cos(phi+!pi/2),sin(phi+!pi/2),replicate(0.0,1,n)]
Btheta_vec = rebin(Btheta,3,n) * [cos(theta+!pi/2)*cos(phi),cos(theta+!pi/2)*sin(phi),sin(theta+!pi/2)]
Bfld       = Bphi_vec + Btheta_vec

; the magnetic axis
Rm = R0+shafr

; finally we convert the normalised minor radius sqrt(psi) to the normalised flux coordinates psi
psi = psi^2

; make the structure
equi = {Bfld:Bfld, psi:psi, theta:theta, Rm:Rm}
; return the result
return, equi

end