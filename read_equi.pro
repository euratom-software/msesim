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

function read_equi, pos, file

; this function reads the data from the equilibrium file and through
; interpolation gets the B-field, normalised flux coordinate, poloidal angle at pos
; and the magnetic axis 
;
; - input:  pos     [3xn] vector-array containing the xyz-coordinates of the points
;           file    the filename of the equilibrium file. This should be an IDL-format file
;                   containing the variables: R, Z, Bfld, Rm and fluxcoord. Such a file can be 
;                   created with the routine 'create_equifile.pro', which simply reads in
;                   the EFIT-equilibrium for a given shot with 'read_flux' and  gets from
;                   that the B-field and the normalised flux coordinates at a given time.
;
; -returns: equi    a structure with fields: Bfld  a [3xn] vector-array with the xyz-components of the B-field at position 'pos'
;                                            psi   a [nx1] array with values in between 0 (magn. axis) and 1 (LCFS) that
;                                                  indicates on which flux surface 'pos' lies
;                                            theta a [nx1] array with values in between 0 and 2*pi that indicates the poloidal angle
;                                            Rm    a scalar that contains the major radius of the magnetic axis

;  v1.0 mdebock, 10/07/2007
;

; read in the file
print,'filename',file
restore, file
print,'equifile:',file

;get the R and Z values of pos
Rpos = sqrt(pos[0,*]^2 + pos[1,*]^2)
Zpos = pos[2,*]
n    = n_elements(Rpos)

; the interpolation indices (Ri,Zi) for these positions are
Ri = (Rpos-min(R))/(max(R)-min(R)) * (n_elements(R)-1)
Zi = (Zpos-min(Z))/(max(Z)-min(Z)) * (n_elements(Z)-1)

; interpolate the B-field components and the flux coordinates
BR   = reform(interpolate(Bfld[0,*,*],Ri,Zi),1,n)
BZ   = reform(interpolate(Bfld[1,*,*],Ri,Zi),1,n)
Bphi = reform(interpolate(Bfld[2,*,*],Ri,Zi),1,n)

psi  = reform(interpolate(transpose(fluxcoord),Ri,Zi),1,n)

; create the B-field vector in xyz-coordinates
Bfld = fltarr(3,n)
phi        = atan(pos[1,*],pos[0,*])					; the toroidal angle of pos
Bphi_vec   = rebin(Bphi,3,n)  * [-sin(phi),cos(phi),replicate(0.0,1,n)]	; the toroidal field vector
Btheta_vec = [BR*cos(phi),BR*sin(phi),BZ]				; the poloidal field vector
Bfld       = Bphi_vec + Btheta_vec					; the total field

; get tje poloidal angle
theta = atan(BR,BZ)-!pi/2
idx   = where(theta lt 0.0, count)
if (count ne 0) then theta(idx) += 2*!pi

; make the structure

equi = {Bfld:Bfld, psi:psi, theta:theta, Rm:Rm}
; return the result

return, equi
end
