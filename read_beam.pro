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

function read_beam, pos, Bv, file

; this function reads the data from the beam file in through
; interpolation gets the beam density, electron density and emission rate at pos 
;
; - input:  pos     [3xn] vector-array containing the xyz-coordinates of the points
;           Bv      [3x1] vector containing the xyz-coordinates of the beam vector
;           file    the filename of the beam file. This should be an IDL-format file
;                   containing the variables: xc, yc, zc, neutral, nel_space and ratebes
;
;  v1.0 mdebock, 07/06/2007
;
;  The interpolation is done to all the gridpoints at once, instead of to just one grid point
;  (so this read-function can be taken out of the loop!)
;  v1.1 mdebock, 25/06/2007
;

; read in the file
print,'filename in read beam',file
restore, file

; The coordinate system is such that the origin is the same as the
; MAST coordinate system and the z-vector is the same as the machine z-vector,
; but the y-vector equals the beam vector. The x-axis follows the right-hand rule.
; Also the unit is cm instead of m

X1 = [Bv[1],-Bv[0],Bv[2]]
; the the coordinates of pos in this coordinate system are:
pos1 = coordtrans(pos,[[0,0,0],[X1]])*100

; the interpolation indices (xi,yi,zi) for pos1 are
xi = (pos1[0,*]-min(xc))/(max(xc)-min(xc)) * (n_elements(xc)-1)
yi = (pos1[1,*]-min(yc))/(max(yc)-min(yc)) * (n_elements(yc)-1)
zi = (pos1[2,*]-min(zc))/(max(zc)-min(zc)) * (n_elements(zc)-1)

neutral2   = neutral[*,*,*,0,0]
nel_space2 = nel_space[*,*,*,0]
ratebes2   = ratebes[*,*,*,0,0]

Bdens    = interpolate(neutral2,xi,yi,zi)
edens    = interpolate(nel_space2,xi,yi,zi)
emission = interpolate(ratebes2,xi,yi,zi)

; Looking at the numbers, it seems as if beam and electron density are given in 
; m^{-3}, whereas emission is given in photons/cm^3/s. To keep it all logical
; we convert the units of emission to photons/m^3/s
emission = 1e6*emission

beam ={Bdens:Bdens, edens:edens, emission:emission}
return, beam

end
