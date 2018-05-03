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

function calc_gp, l, m, efl, C0, Ck, Cl, p0, ds, nslice
; this function calculates x, y and z coordinates of the grid points
; around the position p0.
;
; - input:  l      : 'horizontal' coordinates of the fibres on the focal plane
;           m      : 'vertical'  coordinates of the fibres on the focal plane
;           efl    : effective focal length of the lens
;           C0     : position of the lens in machine coordinates
;           Ck     : k-vector of the lens (optical axis) in machine coordinates
;           Cl     : l-vector of the lens (horizontal) in machine coordinates
;           p0     : central position of the emission
;           ds     : length of one slice 
;           nslice : the number of slices
;
;  v1.0 mdebock, 01/06/2007
;
;  change from single viewing cone to a set of stacked fibres
;  v2.0. mdebock, 14/06/2007
;
;  calculate all the gridpoints at once using an array of relative coordinates (so this
;  calculation can be taken out of the loop!)
;  v2.1 mdebock, 22/06/2007
;
;  Major changes: Coordinates now better represent the closest packing of the fibres.
;  This also means that the input parametes are changed! The 'old' input parameters
;  'fidx' and 'dc' and 'dr' are gone and replaced by 'cols','img_r' and 'nslice'
;  v3.0 mdebock, 13/09/2007
;
;  Major changes: Coordinates now calculated based on the measured/calibrated fibre positions.
;  Also the 'cylindrical' approach is replace by the more real 'cone' shape (i.e. following the
;  emission vectors that converge towards the collection lens - instead of being parallel to
;  one and other).
;  v4.0 mdebock, 10/07/2008
;
;  Minor change: Input coordinates now as separate arrays and not as a structure
;  (this allows different channels to have a different number of grid points)
;  v4.1 mdebock, 10/01/2014

; the distance from the collection lens to the central emission point p0:
D  = norm(C0-p0)
; the sampling along the line-of-sight hence covers the interval:
;  [D-0.5*(nslice-1)*ds,D+0.5*(nslice-1)]

; set up an array to store the gridpoint data
nfib    = n_elements(l)      		; total number of fibres
gp_calc = fltarr(3,nfib*nslice,2)	; gp_calc[*,*,0] contains the gridpoint coordinates,
					; gp_calc[*,*,1] contains the emission vectors

; the emission vectors in the (k,l,m)-coordinate system of the lens are given by
if nfib gt 1. then klm = [replicate(efl,1,nfib), transpose(l), transpose(m)] else klm = [efl,l,m]
klm = klm/rebin(sqrt(klm[0,*]^2 + klm[1,*]^2 + klm[2,*]^2),3,nfib)
; in the (k,l,m)-coordinate system of the lens, the x- and y-axis of the machine coordinate system are:
X1  = coordtrans([1.,0.,0.],[[0.,0.,0.],[Ck],[Cl]])
Y1  = coordtrans([0.,1.,0.],[[0.,0.,0.],[Ck],[Cl]])
; hence the emission vectors in (x,y,z)-machine coordinates are
vec = coordtrans(klm,[[0.,0.,0.],[X1],[Y1]])

; loop through the number of slices and get the gridpoint coordinates
for i=0,nslice-1 do begin
  dist = D-0.5*(nslice-1)*ds + i*ds	; distance gridpoint to the lens for this slice
  gp   = rebin(C0,3,nfib) - dist*vec	; gridpoint coordinates at this distance (distance is
					; subtracted, because the emission vectors point away from the
					; gridpoints and towards the collection lens)
  u = i*nfib
  v = (i+1)*nfib - 1 
  gp_calc[*,u:v,0] = gp			; store gridpoint coordinates
  gp_calc[*,u:v,1] = vec		; and emission vectors
endfor

; return the result
return, gp_calc

end