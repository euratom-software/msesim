pro plot_beam, B0, Bv, Br, Bw
; plots a wireframe cylindre, with a radius of 2* 1/e width of the beam,
; in the currently active window
;
; REMARK: plot_axis (that sets up the window, the rotation, the zoom-factor, ...) 
;         should have been called before calling plot_beam!
;

; settings for the wireframe:
nlines   = 25
ncircles = 12
radius   = Bw
color    = truecolor('darkgray')


; a knot of the wire frame is determined by the point on the beam axis
; and the vector from the beam axis to the 'radius'. In the coordinate system
; that has Bv as X-axis, these vectors are given by:
phi  = transpose(findgen(nlines)*2*!pi/(nlines-1) )
Bwv  = [replicate(0,1,nlines),radius*cos(phi),radius*sin(phi)]
; In this coordinate system the original X-axis [1,0,0] and Y-axis [0,1,0] are:
X1     = coordtrans([1,0,0],[[0,0,0],[Bv]])
Y1     = coordtrans([0,1,0],[[0,0,0],[Bv]])
; Transforming Bwv back into the original coordinate system:
Bwv = coordtrans(Bwv,[[0,0,0],[X1],[Y1]])


; plot the circles around the beam
for k=0,ncircles-1 do begin
  lng = Br[0]+k*(Br[1]-Br[0])/(ncircles-1)
  Bpt = B0+lng*Bv+Bwv[*,0]
  plots, Bpt[0], Bpt[1], Bpt[2], color=color, /T3D, noclip=0
  for l=1,nlines-1 do begin
    Bpt = B0+lng*Bv+Bwv[*,l]
    plots, Bpt[0], Bpt[1], Bpt[2], color=color, /T3D, /continue, noclip=0
   endfor
endfor

; plot the straight lines along the beam
for l=0,nlines-1 do begin
  Bpt = B0+Br[0]*Bv + Bwv[*,l]
  plots, Bpt[0], Bpt[1], Bpt[2], color=color, /T3D, noclip=0
  Bpt = B0+Br[1]*Bv + Bwv[*,l]
  plots, Bpt[0], Bpt[1], Bpt[2], color=color, /T3D, /continue, noclip=0
endfor

end
