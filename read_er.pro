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

function read_Er, pos, Erfile,$
                  equifile,R0,a,shafr,elong,Bphi,q0,qwall,qindex,Bp0,Bpa,Bpindex

; this function reads in the electrostatic Er-field as function of major radius R from a text file (Erfile),
; and converts it to the electrostatic Er-field vector ([3x1] vector) at the given points (pos),
; for the given B-field (equifile, R0,a, etc.)
;
; - input:  pos     [3xn] vector-array containing the xyz-coordinates of the grid points
;           Erfile  filename of the file containing Er as function of R
;
;           equifile file name of the equilibrium (if "none" then calculate the equilibrium)     \
;           R0       major radius of the geometrical centre of the tokamak (also centre LCFS)    |
;           a        minor radius (i.e. radius LCFS)                                             |
;           shafr    Shafranov shift (i.e. magnetic axis at R0+shafr)                            |   needed to convert the R
;           elong    elongation (i.e. height of the plasma is 2*a*elong)                         |>  from the Er-fle into flux
;           Bphi     Toroidal field at R0                                                        |   flux coordinates and to get
;           q0       q at the magnetic axis                                                      |   the radial unit vector
;           qa       q at the plasma edge (LCFS)                                                 |
;           qindex   q = q0 + (qa-q0)*sqrt(psi)^qindex                                           /
;
;
; -returns: Er      [3xn] vector-array containing the Er-field at the each point [Ex,Ey,Ez]
;
;  v1.0 mdebock, 22/04/2008
;


; read in the Er-file
data  = read_txt(Erfile)
data  = reform(data,2,round(0.5*n_elements(data)))
R_in  = data[0,*]
Er_in = data[1,*]
n_in  = n_elements(R_in)

; find the psi-coordinates for R_in
if strcmp('none',equifile,/fold_case) then begin
  equi  = calc_equi([R_in,replicate(0.0,2,n_in)],$          ; if there is no equilibrium file than calculate
                    R0,a,shafr,elong,$                      ; then use the build-in equilibrium model
                    Bphi,q0,qwall,qindex,$
                    Bp0,Bpa,Bpindex)
endif else begin
  equi  = read_equi([R_in,replicate(0.0,2,n_in)], equifile) ; else: get the equilibrium from the file
endelse
psi_in = equi.psi

; sort psi and get rid of duplicate psi
idx    = sort(psi_in)
psi_in = psi_in[idx]
Er_in  = Er_in[idx]
idx    = uniq(psi_in)
psi_in = psi_in[idx]
Er_in  = Er_in[idx]
psi_in = transpose(psi_in)
Er_in  = transpose(Er_in)

; get number of the number of positions
npos = n_elements(pos[0,*])
; get the equilibrium at these positions
if strcmp('none',equifile,/fold_case) then begin
  equi  = calc_equi(pos,$                  ; if there is no equilibrium file than calculate
                    R0,a,shafr,elong,$     ; then use the build-in equilibrium model
                    Bphi,q0,qwall,qindex,$
                    Bp0,Bpa,Bpindex)
endif else begin
  equi  = read_equi(pos, equifile)         ; else: get the equilibrium from the file
endelse
Bfld = equi.Bfld
psi  = equi.psi

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

; Now interpolate the Er-field from psi_in to psi
Er = interpol(Er_in,psi_in, psi)
; and make it a vector
Er = rebin(Er,3,npos)*ur

; and return
return, Er

end