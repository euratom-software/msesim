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

function calc_Er, pos, Bfld, psi, Ermax, Eridx

; this function calculates the electrostatc Er-field ([3x1] vector Er) for the given points (pos),
; the given B-field (Bfld) and the given Er-field parameters (Ermax and Eridx)
;
; - input:  pos     [3xn] vector-array containing the xyz-coordinates of the grid points
;           Bfld    [3xn] vector-array containing the B-field at the each point [Bx,By,Bz]
;           psi     [1xn] array  containing the normalised flux coordinate psi at each point
;           Ermax   Max. in Er-field at sqrt(psi)=0.5
;           Eridx   Er(psi) = Ermax*(1 - (2*abs(0.5-sqrt(psi)))^Eridx)
;
; -returns: Er      [3xn] vector-array containing the Er-field at the each point [Ex,Ey,Ez]
;
;  v1.0 mdebock, 22/04/2008
;

; get number of the number of positions
npos = n_elements(pos[0,*])
; get X, Y, Z and R arrays
X = pos[0,*]
Y = pos[1,*]
Z = pos[2,*]
R = sqrt(pos[0,*]^2 + pos[1,*]^2)
; get Bx,By,Bz arrays
Bx = Bfld[0,*]
By = Bfld[1,*]
Bz = Bfld[2,*]

; rewrite the B-field in terms of Bpol and Bphi:
Bpol = [ X/R^2*(X*Bx + Y*By), Y/R^2*(X*Bx + Y*By), Bz]
Bphi = [-Y/R^2*(X*By - Y*Bx), X/R^2*(X*By - Y*Bx), replicate(0.0,1,npos) ]

; the radial unit vector at each point is given by (Bphi x Bpol)/|Bphi x Bpol|:
ur   = [Bphi[1,*]*Bpol[2,*] - Bphi[2,*]*Bpol[1,*], $
        Bphi[2,*]*Bpol[0,*] - Bphi[0,*]*Bpol[2,*], $
        Bphi[0,*]*Bpol[1,*] - Bphi[1,*]*Bpol[0,*]]
norm = rebin(sqrt(ur[0,*]^2 + ur[1,*]^2 + ur[2,*]^2),3,npos)
ur   = ur/norm
; the above radial unit vector could be pointing inwards instead of outwards
; (if the plasma current and toroidal field are parallel, instead of anti-parallel)
; the way to find out: if Z>0 then ur.ez>0, if Z<0 then ur.ez<0
posZ = where(Z gt 0,count)
if count ne 0 then begin
  if dotp([1.0,0.0,0.0],ur[*,posZ[0]]) lt 0.0 then ur=-ur
endif else begin
  negZ = where(Z lt 0,count)
  if count ne 0 then begin
    if dotp([1.0,0.0,0.0],ur[*,negZ[0]]) gt 0.0 then ur=-ur
  endif
endelse

; Now calculate the Er-field as function of psi
Er = Ermax*(1 - (2*abs(0.5-sqrt(psi)))^Eridx)

Er = reform(Er,1,npos)	; just to make sure that Er is a [1xn] vector
Er = rebin(Er,3,npos)*ur
zeroidx = where(psi gt 1.0, count)
if count ne 0 then Er[*,zeroidx]=replicate(0.0,3,count)

; and return
return, Er

end