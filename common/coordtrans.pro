function coordtrans, V, N 
; function applies a coordinate transformation
; - input:  V  vector-array [3xn] in the old xyz-coordinate system, each row is one point
;           N  can be a [3x2] matrix:
;              the first row gives the origin of the new coordinate system given in the old coordinates.
;              the second row gives the vector describing the new x-axis given in the old coordinates.
;              This new x-axis is 'reached' by rotation around the old z-axis and subsequently around the new
;              y-axis. As a result the new y-axis still lies in the old xy-plane.
;           N  can be a [3x3] matrix:
;              the first row gives the origin of the new coordinate system given in the old coordinates.
;              the second row gives the vector describing the new x-axis given in the old coordinates.
;              the third row gives the vector describing the new y-axis given in the old coordinates. 
;              This new x-axis is 'reached' by rotation around the old z-axis and subsequently around the intermediate
;              y-axis. Finally the new y-axis is reached by rotating around the new x-axis. One should make sure that
;              the vectors for the new x and y axes are perpendicular!
; - returns a vector-array [3xn] containing the xyz-coordinates of the rows of V in the new coordinate system
;
; REMARK: this a 3D transformation. for 2D: set V=[Vx,Vy,0], N = [[Nxx,Nxy,0],[Nyx,Nyy,0],[Nox,Noy,0]] 
; 
;  v1.0 mdebock, 01/06/2007
;
;  v1.1 mdebock, 20/06/2007 : extended the function such that it can handle vector-arrays instead of just 
;                             one vector (and this without using for-loops)
;

; transpose V (such that is becomes an [nx3]-array)
NV = transpose(V)

; get the Nx vector (new x-axis) and the No vector (new origin)
No = N[*,0]
Nx = N[*,1]
if n_elements(N[1,*]) eq 3 then begin
  Ny = N[*,2]
  Ny = Ny/norm(Ny)
endif


; apply the translation
No = transpose(No)			; make No into a [nx3]-array,
No = rebin(No, n_elements(NV[*,0]), 3)	; just like NV
NV = NV - No


; make sure the Nx and Ny vectors are unity vectors
Nx = Nx/norm(Nx)


; the find the rotation angles: first around the old z-axis, then around the intermediate y-axis
angz = atan(Nx[1],Nx[0])
angy = atan(Nx[2],norm(Nx[0:1]))

; the rotation matrices are
A = [[ cos(angz), sin(angz),    0     ],$	; around z
     [-sin(angz), cos(angz),    0     ],$
     [    0     ,    0     ,    1     ]]

B = [[ cos(angy),    0     , sin(angy)],$	; around (intermediate) y
     [    0     ,    1     ,    0     ],$
     [-sin(angy),    0     , cos(angy)]]


if n_elements(N[1,*]) eq 3 then begin
  ; apply the rotations around z and the intermediate on Ny
  Ny = A##Ny
  Ny = B##Ny
  Ny = transpose(Ny) 
  angx = atan(Ny[2],Ny[1])
  ; rotation matrix around the new x
  C = [[    1     ,    0     ,    0     ],$
       [    0     , cos(angx), sin(angx)],$
       [    0     ,-sin(angx), cos(angx)]]
endif else begin
  ; rotation matrix around the new x
  C = [[    1     ,    0     ,    0     ],$
       [    0     ,    1     ,    0     ],$
       [    0     ,    0     ,    1     ]]
endelse

; apply the rotation
NV = A##NV
NV = B##NV
NV = C##NV

; now NV is a [nx3] array, transpose it to a [3xn] vector-array
NV = transpose(NV)

; and return the result
return, NV
end