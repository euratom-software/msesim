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

pro test_Er, calc = calc
; tokamak and field parameters
R0    = 0.83
a     = 0.60
Bphi  = -0.5
q0    = 1.0
qa    = 3.0
qidx  = 3.0
Bp0   = 0.0
Bpa   = 0.0
Bpidx = 1.0
shafr = 0.05
elong = 2.2
Ermax = 2e5
Eridx = 2.0

equifile= 'equi/equi_16246_t0.25.xdr'
;equifile= 'equi/equi_18501_t0.29_TFonly.xdr'

Erfile='Er/Er_test.dat';'none'

; psi and theta
npsi   = 15
ntheta = 40
psi    = findgen(npsi)/(npsi-1)
theta  = findgen(ntheta)*2*!pi/(ntheta-1)
phi    = theta 

; plot a poloidal cross section, with the (poloidal) field indicated
pseudocol,/on
window,0, xsize=1000, ysize=800
!p.background=9
!p.multi = [0,2,2]
plot, [0, 1.05*(R0+a)],[0,0], color=0, thick=2,/isotropic,$
      xs=1,xr=[-0.05,1.1*(R0+a)],ys=1,yr=[-1.1*a*elong, 1.1*a*elong]
xyouts, 1.05*(R0+a), -0.1*(a*elong), 'R', color=0,charsize=1.5
oplot, [0,0], [-1.05*a*elong, 1.05*a*elong], color=0, thick=2
xyouts, 0.05*(R0+a), 1.0*(a*elong), 'Z', color=0,charsize=1.5

pos = fltarr(3,npsi*ntheta)
m   = 0
for k=0,npsi-1 do begin
for l=0,ntheta-1 do begin
  R        = R0 + shafr*(1 - psi[k]^2) + psi[k]*a*cos(theta[l])
  Z        = psi[k]*a*elong*sin(theta[l])
  pos[*,m] = [R,0,Z]
  m++
endfor
endfor

; get equilibrium
if keyword_set(calc) then begin
  equi = calc_equi(pos, R0, a, shafr, elong, Bphi, q0, qa, qidx,Bp0,Bpa,Bpidx)
endif else begin
  equi = read_equi(pos, equifile)
endelse


Bfld   = equi.Bfld
psi2   = equi.psi
theta2 = equi.theta

; calculate Er
if Erfile eq 'none' then begin
  Er = calc_Er(pos, Bfld, psi2, Ermax, Eridx)
endif else begin
  Er = read_Er(pos, Erfile,equifile,R0,a,shafr,elong,Bphi,q0,qa,qidx,Bp0,Bpa,Bpidx)
endelse
; and plot it
for m=0,npsi*ntheta-1 do begin
  R  = pos[0,m]
  Z  = pos[2,m]
  BR = Bfld[0,m]
  BZ = Bfld[2,m]
  arrow, R,Z, R+BR,Z+BZ, /data, color=0,thick=1, hthick=1, hsize=(!D.X_SIZE/200.),/solid
  arrow, R,Z, R+Er[0,m]*2e-6,Z+Er[2,m]*2e-6, /data, color=1,thick=1, hthick=1, hsize=-0.2,/solid
end

; plot a top view, with the (toroidal) field indicated
plot, [-1.05*(R0+a), 1.05*(R0+a)],[0,0], color=0, thick=2,/isotropic,$
      xs=1,xr=[-1.1*(R0+a),1.1*(R0+a)],ys=1,yr=[-1.1*(R0+a), 1.1*(R0+a)]
xyouts, 1.05*(R0+a), -0.1*(R0+a), 'X', color=0,charsize=1.5
oplot, [0,0], [-1.05*(R0+a), 1.05*(R0+a)], color=0, thick=2
xyouts, 0.05*(R0+a), 1.0*(R0+a), 'Y', color=0,charsize=1.5

pos = fltarr(3,2*npsi*ntheta)
m=0
for k=0,npsi-1 do begin
for l=0,ntheta-1 do begin
  R        = R0 + shafr*(1 - psi[k]^2) + psi[k]*a
  X        = R*cos(theta[l])
  Y        = R*sin(theta[l])
  pos[*,m] = [X,Y,0]
  m++
  R        = R0 + shafr*(1 - psi[k]^2) - psi[k]*a
  X        = R*cos(theta[l])
  Y        = R*sin(theta[l])
  pos[*,m] = [X,Y,0]
  m++
endfor
endfor
if keyword_set(calc) then begin
  equi = calc_equi(pos, R0, a, shafr, elong, Bphi, q0, qa, qidx,Bp0,Bpa,Bpidx)
endif else begin
  equi = read_equi(pos, equifile)
endelse

Bfld = equi.Bfld
psi2 = equi.psi
Rm   = equi.Rm
; calculate Er
if Erfile eq 'none' then begin
  Er = calc_Er(pos, Bfld, psi2, Ermax, Eridx)
endif else begin
  Er = read_Er(pos, Erfile,equifile,R0,a,shafr,elong,Bphi,q0,qa,qidx,Bp0,Bpa,Bpidx)
endelse
; and plot it
for m=0,2*npsi*ntheta-1 do begin
  X  = pos[0,m]
  Y  = pos[1,m]
  BX = Bfld[0,m]
  BY = Bfld[1,m]
  arrow,X,Y,  X+BX,Y+BY, /data, color=0,thick=1, hthick=1, hsize=(!D.X_SIZE/200.),/solid
  arrow,X,Y,  X+Er[0,m]*2e-6,Y+Er[1,m]*2e-6, /data, color=1,thick=1, hthick=1, hsize=-0.2,/solid
end

; plot norm(Er) as function of psi2
X = pos[0,*]
Y = pos[1,*]
Ex = Er[0,*]
Ey = Er[1,*]
Ernorm = sqrt(Er[0,*]^2 + Er[1,*]^2 + Er[2,*]^2)
idx = where(sqrt(pos[0,*]^2 + pos[1,*]^2) ge Rm,count) ; plot only for R>Rm
if count ne 0 then begin
  Ernorm = Ernorm[idx]
  psi3   = psi2[idx]
  X      = X[idx]
  Y      = Y[idx]
  Ex     = Ex[idx]
  Ey     = Ey[idx]
endif
Eproj  = X*Ex + Y*Ey		; projection on R-unit vector: positive => outward Er, negative => inward Er
idx = where(Eproj lt 0.0, count)
if count ne 0 then Ernorm[idx] = -Ernorm[idx]
plot,sqrt(psi3),Ernorm,$
     color=0, xtitle='sqrt(psi)', ytitle='Er'

!p.multi=0




end
