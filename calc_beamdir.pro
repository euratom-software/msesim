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

function calc_beamdir, pos, B0, Bv, hf, vf, div, pini

; calculates the unit vectors indicating the velocity directions resulting from
; each pini-hole and their propability.
;
; - input:  pos      [3xn] vector-array containing the xyz-coordinates of the gridpoints
;           B0       coordinates of the beam source/pini (NOT THE BEAM DUCT)
;           Bv       the beam vector
;           hf       distance between the beam source and the horizontal focus
;           vf       distance between the beam source and the vertical focus
;           div      half 1/e width of the divergence
;           pini     [2xm] array with the yz-coordinates of the pini-holes in the beam
;                    coordinate system. 
;
; -returns: vdir     [nx1]-structure array with following fields:
;                       - avvel: [3x1] average (unit) velocity vector
;                       - vel  : [3xm] vector-array (unit) velocity direction vectors from each pini-hole
;                       - p    : [1xm] array with the propability for each velocity direction
;
;  v1.0 mdebock, 19/07/2007
;
;  v1.1 mdebock, 25/07/2007 : hacky way to include the old beam model, in order to find out whether
;                             the differences are real
;

; this is the hacky, hard-coded way to include the old beam model and the plot the velocity distribution
oldbeam=0
an_div =4
seg_div=5
doplot =0


if ~oldbeam then begin
  ; number of grid points and the number of pini-holes
  npos  = n_elements(pos[0,*])
  npini = n_elements(pini[0,*])

  ; in the coordinate system that has the beam-vector Bv as X-axis and
  ; B0 as origin the coordinates of the gridpoints are:
  pos1 = coordtrans(pos,[[B0],[Bv]])
  ; the original X- and Y-axis become:
  X1 = coordtrans([1,0,0],[[0,0,0],[Bv]])
  Y1 = coordtrans([0,1,0],[[0,0,0],[Bv]])

  ; get the unit direction vector for each pini-hole
  vert = fltarr(3,npini)
  horz = fltarr(3,npini)
  ; first vertical:
  vert = rebin([vf,0.0,0.0],3,npini) - [fltarr(1,npini),fltarr(1,npini),pini[1,*]]
  ; then horizontal:
  horz = rebin([hf,0.0,0.0],3,npini) - [fltarr(1,npini),pini[0,*],fltarr(1,npini)]
  ; combine vertical and horizontal
  vpini      = fltarr(3,npini)
  vpini[0,*] = vert[0,*]
  vpini[1,*] = vert[0,*]/horz[0,*] * horz[1,*]
  vpini[2,*] = vert[2,*]
  ; and make it unit vectors
  vpini = vpini / rebin(sqrt(vpini[0,*]^2 + vpini[1,*]^2 + vpini[2,*]^2),3,npini)

  ; create structure array
  struct    = {avvel:[0.0,0.0,0.0], vel:fltarr(3,npini), p:fltarr(1,npini)}
  structarr = replicate(struct,npos)
  ; now loop through the grid points
  for k=0l,npos-1 do begin
    ; get the unit vectors from each pini-hole to the gridpoint: i.e. the velocity directions 
    vel1 = rebin(pos1[*,k],3,npini) - [fltarr(1,npini),pini[0,*],pini[1,*]]
    vel1 = vel1 / rebin(sqrt(vel1[0,*]^2 + vel1[1,*]^2 + vel1[2,*]^2),3,npini)

    ; get the angle xi between the pini direction and vel
    dp  = vpini[0,*]*vel1[0,*] + vpini[1,*]*vel1[1,*] + vpini[2,*]*vel1[2,*]
    idx = where(dp gt 1.0,count)		; can happen due to round-off error
    if count ne 0 then dp[idx]=1.0
    idx = where(dp lt -1.0,count)	; can happen due to round-off error
    if count ne 0 then dp[idx]=-1.0
    xi  = acos(dp)
    ;print,xi[0]
    ; get the propability
    p   = exp(-xi^2/div^2)
    p   = p/total(p)

    ; get the average velocity direction
    if npini gt 1 then avvel1 = total(rebin(p,3,npini)*vel1,2) else avvel1 = vel1
    avvel1 = avvel1/norm(avvel1)

    ; transform the vectors back into the machine coordinate system
    vel   = coordtrans(vel1,[[0,0,0],[X1],[Y1]])
    avvel = coordtrans(avvel1,[[0,0,0],[X1],[Y1]])

    ; put all data in the structure
    structarr[k].avvel = avvel
    structarr[k].vel   = vel
    structarr[k].p     = p

  endfor

endif else begin

  ; number of grid points
  npos  = n_elements(pos[0,*])

  ; in the coordinate system that has the beam-vector Bv as X-axis and
  ; B0 as origin the coordinates of the gridpoints are:
  pos1 = coordtrans(pos,[[B0],[Bv]])
  ; the original X- and Y-axis become:
  X1 = coordtrans([1,0,0],[[0,0,0],[Bv]])
  Y1 = coordtrans([0,1,0],[[0,0,0],[Bv]])

  ; intialise vector-arrays for both a vertical and horizontal direction
  vert = fltarr(3,npos)
  horz = fltarr(3,npos)

  ; first vertical:
  ; the beam-particle direction depends on whether we are converging to are diverging from the focus
  idx = where(pos1[0,*] ge vf, count)
  if (count ne 0) then begin
    vert[*,idx] = [pos1[0,idx], fltarr(1,count), pos1[2,idx]] - rebin([vf,0.0,0.0],3,count)
  endif
  idx = where(pos1[0,*] lt vf, count)
  if (count ne 0) then begin
    vert[*,idx] = rebin([vf,0.0,0.0],3,count) - [pos1[0,idx], fltarr(1,count), pos1[2,idx]]
  endif

  ; then horizontal:
  ; the beam-particle direction depends on whether we are converging to are diverging from the focus
  idx = where(pos1[0,*] ge hf, count)
  if (count ne 0) then begin
    horz[*,idx] = [pos1[0,idx], pos1[1,idx], fltarr(1,count)] - rebin([hf,0.0,0.0],3,count)
  endif
  idx = where(pos1[0,*] lt hf, count)
  if (count ne 0) then begin
    horz[*,idx] = rebin([hf,0.0,0.0],3,count) - [pos1[0,idx], pos1[1,idx], fltarr(1,count)]
  endif

  ; combine vertical and horizontal
  avvel1      = fltarr(3,npos)
  avvel1[0,*] = vert[0,*]
  avvel1[1,*] = vert[0,*]/horz[0,*] * horz[1,*]
  avvel1[2,*] = vert[2,*]
  ; and make it unit vectors
  avvel1 = avvel1 / rebin(sqrt(avvel1[0,*]^2 + avvel1[1,*]^2 + avvel1[2,*]^2),3,npos)

  ; put it back in the original coordinate system
  avvel = coordtrans(avvel1,[[0,0,0],[X1],[Y1]])


  ; Calculate the divergence directions and the propability
  ; The radius of the 'divergence disk' corresponds with 2 times the 1/e width
  div_r   = sin(2.*div)
  ; The radius corresponding with the 1/e width:
  div_w   = sin(div)
  ; We use the function calc_dp to calculate the coordinates of a number of
  ; points homogeneously distributed over the 'divergence disk'
  div2 = calc_dp(div_r,an_div,seg_div)
  ; the number of points:
  ndiv = n_elements(div2[0,*])
  ; For the propability we use a normaised Gaussian distribution:
  p = exp(-(div2[0,*]^2 + div2[1,*]^2)/div_w^2)
  p = p/total(p)


  ; create structure array
  struct    = {avvel:[0.0,0.0,0.0], vel:fltarr(3,ndiv), p:fltarr(1,ndiv)}
  structarr = replicate(struct,npos)
  ; now loop through the grid points
  for k=0,npos-1 do begin
    ; In the coordinate system with avvel as the x-axis the orignal x- and y-axis are given by:
    X2  = coordtrans([1,0,0],[[0,0,0],[avvel[*,k]]])
    Y2  = coordtrans([0,1,0],[[0,0,0],[avvel[*,k]]])
    ; With X2 and Y2 known we can transform the coordinates of the divergence vectors in the
    ; plane perpendicular to the velocity direction:
    div_vec = coordtrans([replicate(0.0,1,ndiv),div2],[[0.0,0.0,0.0],[X2],[Y2]])

    ; Loop through the beam divergence
    vel = fltarr(3,ndiv)
    for i=0,ndiv-1 do begin
      ; the diverged velocity direction is the normalised sum of the average direction and the divergence vector
      vel[*,i] = avvel[*,k] + div_vec[*,i]
      vel[*,i] = vel[*,i]/norm(vel[*,i])
    endfor

    ; put all data in the structure
    structarr[k].avvel = avvel[*,k]
    ;print, 'average velocity distribution', avvel[*,k]
    structarr[k].vel   = vel
    structarr[k].p     = p
  endfor
endelse


if doplot then begin
  pseudocol
  !p.background=9
  !p.multi = [0,3,3]
  window,0,xsize=850,ysize=850
  for k=0,npos-1 do begin
    plot,[-0.025,0.025],[-0.045,0.045],/nodata, color=0,/isotropic, xs=1, ys=1
    p     = structarr[k].p
    vel   = structarr[k].vel
    avvel = structarr[k].avvel
    color = 10+round(245 * p/max(p))
    plots, vel[1,*], vel[2,*], psym=1, color=color, thick=2

    plots, avvel[1], avvel[2], psym=1, color=0, symsize=3,thick=2;

    oplot,[-0.04,0.04],[0.0,0.0], color=0, linestyle=2
    oplot,[0.0,0.0],[-0.10,0.10], color=0, linestyle=2
  endfor
endif


; finally return the result
return, structarr


end
