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

function calc_dp, rtot, rings, segments, doplot=doplot, noline=noline, anaz=anaz
;
; function calculates a homogeneous distribution of points on a disk. The area
; around each disk-point is equal
;
; - input:  rtot     : the radius of the disk
;           rings    : sets number of rings in which the disk is divided
;           segments : the increment number of segments (ring0:  1, ring1 = 1 + segment, ring2= 1+2*segment, ... )
;
;           the keyword 'anaz' will treat rings as annuli and segments as azimuths thus creating
;           diskpoints in the 'old' way
; 
; -returns: dp   : the xy-coordinates of the disk-points [2xn]
;
;  v1.0 mdebock, 29/06/2007




if keyword_set(noline) then color=0 else color=9

; check whether the input parameters are correct
rings    = round(rings)
segments = round(segments)

if (rtot lt 0) then begin
  print, 'ERROR: The radius cannot be negative!'
  return, -1
endif

if (rings lt 1) then begin
  print, 'ERROR: The number of rings has to be larger than zero!'
  return, -1
endif

if (segments lt 1) then begin
  print, 'ERROR: The number of segments has to be larger than zero!'
  return, -1
endif


; if r=0 => only 1 central point 
if abs(rtot) lt 1e-6 then return, dp = [0.0,0.0]


if not(keyword_set(anaz)) then begin

; array of radii
r    = fltarr(rings)
; array of diskpoints
ntot = 1+total(indgen(rings)*segments)
dp   = fltarr(2,ntot)

; find the central point and the
; radius of the inner circle
r[0] = rtot/sqrt(ntot)
dp[*,0] = [0,0]


; plot if wanted
if keyword_set(doplot) then begin
  pseudocol
  !p.background=0
  phi = findgen(100) * 2*!pi/99.
  plot, r[0]*cos(phi),r[0]*sin(phi), color=color, /isotropic,$
        xs=1, xr=[-1.05*rtot,1.05*rtot],ys=1, yr=[-1.05*rtot,1.05*rtot]
  plots, dp[0,0],dp[1,0],psym=1,symsize=2,thick=2,color=1
endif


; find the other points
j=1
for l=1,rings-1 do begin
  ; radius of the ring is determined by the condition of equal area
  r[l]   = sqrt(segments*l*r[0]^2 + r[l-1]^2)

  ; plot the ring if wanted
  if keyword_set(doplot) then begin
    oplot, r[l]*cos(phi),r[l]*sin(phi), color=color
  endif

  ; the opening angle of the segments on the ring
  dtheta = 2*!pi/(segments*l)
  for m=0,segments*l-1 do begin
    ; determine the central point of the segment
    theta = (0.5 + m)*dtheta
    rdp   = (r[l]+r[l-1])/2
    dp[*,j] = [rdp*cos(theta),rdp*sin(theta)]

    ; plot the segment if wanted
    if keyword_set(doplot) then begin
      plots, dp[0,j],dp[1,j],psym=1,symsize=2,thick=2,color=1
      oplot, [r[l-1]*cos(theta-0.5*dtheta),r[l]*cos(theta-0.5*dtheta)],$
             [r[l-1]*sin(theta-0.5*dtheta),r[l]*sin(theta-0.5*dtheta)],$
             color=color
    endif
    j++
  endfor
endfor

; return the result
return, dp

endif else begin

; array of radii
r    = fltarr(rings+1)
; array of diskpoints
ntot = rings*segments
dp   = fltarr(2,ntot)

; plot if wanted
if keyword_set(doplot) then begin
  pseudocol
  !p.background=0
  phi = findgen(100) * 2*!pi/99.
  plot, r[0]*cos(phi),r[0]*sin(phi), color=color, /isotropic,$
        xs=1, xr=[-1.05*rtot,1.05*rtot],ys=1, yr=[-1.05*rtot,1.05*rtot], /nodata
  plots, 0,0,psym=1,symsize=2,thick=2,color=color
endif


; find the other points
j=0
r[0] = 0.0
dtheta = 2*!pi/segments
for l=1,rings do begin
  ; radius of the ring is determined by the condition of equal area
  rdp = sqrt(2*l-1)*rtot/sqrt((2*rings))
  r[l] = sqrt(2*l)*rtot/sqrt((2*rings))

  ; plot the ring if wanted
  if keyword_set(doplot) then begin
    oplot, r[l]*cos(phi),r[l]*sin(phi), color=color
  endif

  ; the opening angle of the segments on the ring
  for m=0,segments-1 do begin
    ; determine the central point of the segment
    theta = m*dtheta
    dp[*,j] = [rdp*cos(theta),rdp*sin(theta)]

    ; plot the segment if wanted
    if keyword_set(doplot) then begin
      plots, dp[0,j],dp[1,j],psym=1,symsize=2,thick=2,color=1
      oplot, [r[l-1]*cos(theta-0.5*dtheta),r[l]*cos(theta-0.5*dtheta)],$
             [r[l-1]*sin(theta-0.5*dtheta),r[l]*sin(theta-0.5*dtheta)],$
             color=color
    endif
    j++
  endfor
endfor

; return the result
return, dp

endelse

end