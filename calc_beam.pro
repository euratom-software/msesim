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

function calc_beam, pos,B0,Bv,chi,w0,Bdens0,edens0,Qion,Qemit,$
                    equifile,R0,a,shafr,elong,Bphi,q0,qwall,qindex,Bp0,Bpa,Bpindex

; this function calculates the beam density, the electron density and the emission rate at positions pos
;
; - input:  pos     [3xn] vector-array containing the xyz-coordinates of the grid points
;           B0      coordinates of the beam port
;           Bv      the beam vector
;           chi     the beam divergence
;           w0      the beam width (half 1/e) at the beam port
;           Bdens0  beam density outside the plasma (at R=R0+a) 
;           edens0  electron density at the magetic axis
;           Qion    ionisation rate of beam particles
;           Qemit   emission rate of beam particles
;
;           equifile file name of the equilibrium (if "none" then calculate the equilibrium)     \
;           R0       major radius of the geometrical centre of the tokamak (also centre LCFS)    |
;           a        minor radius (i.e. radius LCFS)                                             |
;           shafr    Shafranov shift (i.e. magnetic axis at R0+shafr)                            |
;           elong    elongation (i.e. height of the plasma is 2*a*elong)                         |>  needed for [x,y,z] to
;           Bphi     Toroidal field at R0                                                        |   flux coordinates conversion
;           q0       q at the magnetic axis                                                      |
;           qa       q at the plasma edge (LCFS)                                                 |
;           qindex   q = q0 + (qa-q0)*sqrt(psi)^qindex                                           /
;
; -returns: beam    a structure with fields: Bdens    beam density at pos      [1xn]-array
;                                            edens    electron density at pos  [1xn]-array
;                                            emission emission rate = Bdens*edens*Qemit   [1xn]-array
;
;  v1.0 mdebock, 06/06/2007
;
;  calculate the beam at all the gridpoints at once, instead of calculating it for just one grid point
;  (so this calculation can be taken out of the loop!)
;  v1.1 mdebock, 22/06/2007
;
;  inclusion of equifile
;  v1.2 mdebock, 22/04/2008


; the beam model is a simple one:
; - we assume a parabolic density profile with edens=edens0 for psi=0 and edens=0 for psi=1
; - it is assumed beam particles are only lost on the beam axis: 
;             d Bdens  = - edens * Bdens * Qion * dl
;
;        <=>  Bdens(l) = Bdens0 * exp[ - Qion * integral(0->l)(edens(l)dl) ]
;
;   where l is the distance from the entry point (where beam axis crosses the LCFS
;
; - perpendicular to the beam axis the profile is Gaussian, with a 1/e width:
;             w(k) = w0 + k*sin(chi)
;   where k is the distance to the beam port B0
;
; - the emission rate is taken to be 2e13 (in the same order as what M. Turnyanskiy predicts)


; Numerical integration setting:
nl = 100    ; number of numerical integration points

; number of pos:
npos  = n_elements(pos[0,*])
; initialise the arrays that will hold the final results
Bdens    = fltarr(1,npos)
Bdensc   = fltarr(1,npos)	; the central beam density 
edens    = fltarr(1,npos)
emission = fltarr(1,npos)


; first thing to do is to find out the distance 'd' from pos to the neutral beam axis
; and the distance 'k' from B0 to the 'footprint' of pos on the neutral beam axis
B0pos_v = pos-rebin(B0,3,npos)		; vectors from B0 to pos
d       = sqrt( (Bv[1]*B0pos_v[2,*]-Bv[2]*B0pos_v[1,*])^2$
               +(Bv[0]*B0pos_v[2,*]-Bv[2]*B0pos_v[0,*])^2$
               +(Bv[0]*B0pos_v[1,*]-Bv[1]*B0pos_v[0,*])^2)
n       = sqrt(B0pos_v[0,*]^2 + B0pos_v[1,*]^2 + B0pos_v[2,*]^2)
k       = fltarr(1,npos)
idx     = where(n ge d, count)		; due to floating point rounding errors it is possible that n<d,
if (count ne 0) then begin		; of course this is nonsense and just means that n=d and k=0
  k[idx] = sqrt(n[idx]^2 - d[idx]^2)	; Pythagoras to find distance k along the beam axis from B0 to the footprint
endif

; we now try to find the distance m from B0 to the entry point. I.e. where R=R0+a 
; <=> (B0[0] + m*Bv[0])^2 + (B0[0] + m*Bv[0])^2 = (R0+a)^2
; <=> (Bv[0]^2 + Bv[1]^2)*m^2 + 2*(B0[0]*Bv[0] + B0[1]*Bv[1])*m  + [(B0[0]^2 + B0[1]^2  - (R0+a)^2] = 0
coeff = [B0[0]^2 + B0[1]^2  - (R0+a)^2,$
         2*(B0[0]*Bv[0] + B0[1]*Bv[1]),$
         Bv[0]^2 + Bv[1]^2             ]
m     = real_part(fz_roots(coeff))
m     = min(m)				; because it's the entry point we're looking for, not the exit point
; the distance l between the entry point and the footprint then is:
l     = k-m

; if there would be no plasma, then the central beam density would be everywhere the same:
Bdensc[*] = Bdens0

; when there is a plasma, we need the integrate over the electron density along the
; beam axis. So we should find the largest 'l' (i.e. the point deepest in the plasma)
; and find the electron density at 'nl' points in between 0 and 'l'
idx = where(l gt 0.0, count)
if (count ne 0) then begin
  ; all the positive l's and the maximum l
  lpos = l[idx]
  lmx = max(lpos)

  ; step in l
  dl    = lmx/nl
  lbin  = (0.5 + findgen(1,nl))*dl

  ; the xyz-coordinates of the points on the beam
  Bpt = rebin(B0,3,nl) + (m + rebin(lbin,3,nl)) * rebin(Bv,3,nl)

  ; find the psi-coordinates
  if strcmp('none',equifile,/fold_case) then begin
  equi  = calc_equi(Bpt,$			; if there is no equilibrium file than calculate
                    R0,a,shafr,elong,$		; then use the build-in equilibrium model
                    Bphi,q0,qwall,qindex,$
                    Bp0,Bpa,Bpindex)
  endif else begin
    equi  = read_equi(Bpt, equifile)		; else: get the equilibrium from the file
  endelse
  psi   = equi.psi

  ; and get the integrated electron density at these points (parabolic density profile assumed)
  intdens = total(edens0*(1-sqrt(psi)^2),/cumulative)*dl

  ; the l-bin in which l falls.
  lidx = ceil((lpos-dl)/dl)

  ; we now know the beam density at the footprint
  Bdensc[idx] =  Bdens0*exp(-Qion*intdens[lidx])
endif

; next thing is to find out the 1/e width of the beam at the footprint.
; that depends on the width at B0 and the beam divergence:
w = w0+ k*sin(chi)
; the beam density at pos then is:
Bdens = (Bdensc*w0)/(w) * exp(-d^2/w^2)

; the electron density at pos:
; calculate the equilibrium to find the psi
;equi  = calc_equi(pos,R0,a,shafr,elong,Bphi,q0,qwall,qindex,Bp0,Bpa,Bpindex)

if strcmp('none',equifile,/fold_case) then begin
equi = calc_equi(pos,$
                 R0,a,sharf,elong,$
                 Bphi, q0, qwall,qindex,$
                 Bp0, Bpa, Bpindex)
end if else begin

equi = read_equi(pos, equifile)
psi   = equi.psi

; get the electron density (parabolic density profile assumed)
edens = edens0*(1-psi^2)

; the last thing is the emission
; REMARK: to avoid floating point overflows, we convert Bdens and edens to 1e12 m^(-3), 
;         and Qemit to 1e-20 photons*m^3/s/sr
Bdenstmp = Bdens*1e-12
edenstmp = edens*1e-12
Qemittmp = Qemit*1e20
emission = Bdenstmp * edenstmp * Qemittmp

; make the structure and return it
beam ={Bdens:Bdens, edens:edens, emission:emission}
return, beam

end
