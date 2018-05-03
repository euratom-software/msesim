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

pro plot_axis, centre, zoom, rotation, ttl, lth

; plots the 3D axis
; the axis should be plotted before gridpoints, neutral beam, magnetic field ... are plotted!


; define the plot area (the drawn axes are 5% shorter at the plot area) 
xr=[centre[0]-1.0/zoom, centre[0]+1.0/zoom]
yr=[centre[1]-1.0/zoom, centre[1]+1.0/zoom]
zr=[centre[2]-1.0/zoom, centre[2]+1.0/zoom]

; erase the plot (i.e. set background colour)
erase

; set character size if not set
oldcharsize = !p.charsize
if !p.charsize eq 0 then !p.charsize=1

; set up the plot ranges and rotation (viewing angle)az = rotation[0] mod 360.
az = rotation[0] mod 360.
ax = abs(rotation[1] mod 90.)
if az lt 0 then az += 360.
scale3, xr=xr, yr=yr, zr=zr, az=az,ax=ax			

; plot the axes (where to plot them depends on the rotation)
if (az ge 0.) && (az lt 90.) then begin
  plot, xr, yr, /T3D, /nodata, color=truecolor('black'), xtitle='', ytitle='',$
        xs=1, ys=1, xtickname=replicate(' ',30), ytickname=replicate(' ',30), /noerase
  axis, xaxis=0, /T3D, xtitle='X [m]', color=truecolor('black'), xs=1, xr=xr
  axis, yaxis=0, /T3D, ytitle='Y [m]', color=truecolor('black'), ys=1, yr=yr
  axis, zaxis=2, /T3D, ztitle='Z [m]', color=truecolor('black'), zs=1, zr=zr
  axis, zaxis=0, /T3D, ztitle='', color=truecolor('black'), ztickname=replicate(' ',30), zs=1, zr=zr
  axis, zaxis=3, /T3D, ztitle='', color=truecolor('black'), ztickname=replicate(' ',30), zs=1, zr=zr
  axis, xr[1], yr[1], zr[1], /xaxis, /T3D, xtitle='', color=truecolor('black'), xtickname=replicate(' ',30), xs=1, xr=xr
  axis, xr[1], yr[1], zr[1], /yaxis, /T3D, xtitle='', color=truecolor('black'), ytickname=replicate(' ',30), ys=1, yr=yr
  if (zr[0] lt 0) && (zr[1] gt 0) then begin
    plots, xr, [yr[1],yr[1]], [0,0], linestyle=2, color=truecolor('black'), /T3D
    plots, [xr[1],xr[1]], yr, [0,0], linestyle=2, color=truecolor('black'), /T3D
  endif
  if (yr[0] lt 0) && (yr[1] gt 0) then plots, xr, [0,0], [zr[0],zr[0]], linestyle=2, color=truecolor('black'), /T3D
  if (xr[0] lt 0) && (xr[1] gt 0) then plots, [0,0], yr, [zr[0],zr[0]], linestyle=2, color=truecolor('black'), /T3D
endif 

if (az ge 90.) && (az lt 180.) then begin
  plot, xr, yr, /T3D, /nodata, color=truecolor('black'), xtitle='', ytitle='',$
        xs=1, ys=1, xtickname=replicate(' ',30), ytickname=replicate(' ',30), /noerase
  axis, xaxis=1, /T3D, xtitle='X [m]', color=truecolor('black'), xs=1, xr=xr
  axis, yaxis=0, /T3D, ytitle='Y [m]', color=truecolor('black'), ys=1, yr=yr
  axis, zaxis=1, /T3D, ztitle='Z [m]', color=truecolor('black'), zs=1, zr=zr
  axis, zaxis=0, /T3D, ztitle='', color=truecolor('black'), ztickname=replicate(' ',30), zs=1, zr=zr
  axis, zaxis=3, /T3D, ztitle='', color=truecolor('black'), ztickname=replicate(' ',30), zs=1, zr=zr
  axis, xr[1], yr[0], zr[1], /xaxis, /T3D, xtitle='', color=truecolor('black'), xtickname=replicate(' ',30), xs=1, xr=xr
  axis, xr[1], yr[1], zr[1], /yaxis, /T3D, xtitle='', color=truecolor('black'), ytickname=replicate(' ',30), ys=1, yr=yr
  if (zr[0] lt 0) && (zr[1] gt 0) then begin
    plots, xr, [yr[0],yr[0]], [0,0], linestyle=2, color=truecolor('black'), /T3D
    plots, [xr[1],xr[1]], yr, [0,0], linestyle=2, color=truecolor('black'), /T3D
  endif
  if (yr[0] lt 0) && (yr[1] gt 0) then plots, xr, [0,0], [zr[0],zr[0]], linestyle=2, color=truecolor('black'), /T3D
  if (xr[0] lt 0) && (xr[1] gt 0) then plots, [0,0], yr, [zr[0],zr[0]], linestyle=2, color=truecolor('black'), /T3D
endif 

if (az ge 180.) && (az lt 270.) then begin
  plot, xr, yr, /T3D, /nodata, color=truecolor('black'), xtitle='', ytitle='',$
        xs=1, ys=1, xtickname=replicate(' ',30), ytickname=replicate(' ',30), /noerase
  axis, xaxis=1, /T3D, xtitle='X [m]', color=truecolor('black'), xs=1, xr=xr
  axis, yaxis=1, /T3D, ytitle='Y [m]', color=truecolor('black'), ys=1, yr=yr
  axis, zaxis=0, /T3D, ztitle='Z [m]', color=truecolor('black'), zs=1, zr=zr
  axis, zaxis=1, /T3D, ztitle='', color=truecolor('black'), ztickname=replicate(' ',30), zs=1, zr=zr
  axis, zaxis=2, /T3D, ztitle='', color=truecolor('black'), ztickname=replicate(' ',30), zs=1, zr=zr
  axis, xr[0], yr[0], zr[1], /xaxis, /T3D, xtitle='', color=truecolor('black'), xtickname=replicate(' ',30), xs=1, xr=xr
  axis, xr[0], yr[1], zr[1], /yaxis, /T3D, xtitle='', color=truecolor('black'), ytickname=replicate(' ',30), ys=1, yr=yr
  if (zr[0] lt 0) && (zr[1] gt 0) then begin
    plots, xr, [yr[0],yr[0]], [0,0], linestyle=2, color=truecolor('black'), /T3D
    plots, [xr[0],xr[0]], yr, [0,0], linestyle=2, color=truecolor('black'), /T3D
  endif
  if (yr[0] lt 0) && (yr[1] gt 0) then plots, xr, [0,0], [zr[0],zr[0]], linestyle=2, color=truecolor('black'), /T3D
  if (xr[0] lt 0) && (xr[1] gt 0) then plots, [0,0], yr, [zr[0],zr[0]], linestyle=2, color=truecolor('black'), /T3D
endif 

if (az ge 270.) && (az lt 360.) then begin
  plot, xr, yr, /T3D, /nodata, color=truecolor('black'), xtitle='', ytitle='',$
        xs=1, ys=1, xtickname=replicate(' ',30), ytickname=replicate(' ',30), /noerase
  axis, xaxis=0, /T3D, xtitle='X [m]', color=truecolor('black'), xs=1, xr=xr
  axis, yaxis=1, /T3D, ytitle='Y [m]', color=truecolor('black'), ys=1, yr=yr
  axis, zaxis=1, /T3D, ztitle='Z [m]', color=truecolor('black'), zs=1, zr=zr
  axis, zaxis=2, /T3D, ztitle='', color=truecolor('black'), ztickname=replicate(' ',30), zs=1, zr=zr
  axis, zaxis=3, /T3D, ztitle='', color=truecolor('black'), ztickname=replicate(' ',30), zs=1, zr=zr
  axis, xr[0], yr[1], zr[1], /xaxis, /T3D, xtitle='', color=truecolor('black'), xtickname=replicate(' ',30), xs=1, xr=xr
  axis, xr[0], yr[1], zr[1], /yaxis, /T3D, xtitle='', color=truecolor('black'), ytickname=replicate(' ',30), ys=1, yr=yr
  if (zr[0] lt 0) && (zr[1] gt 0) then begin
    plots, xr, [yr[1],yr[1]], [0,0], linestyle=2, color=truecolor('black'), /T3D
    plots, [xr[0],xr[0]], yr, [0,0], linestyle=2, color=truecolor('black'), /T3D
  endif
  if (yr[0] lt 0) && (yr[1] gt 0) then plots, xr, [0,0], [zr[0],zr[0]], linestyle=2, color=truecolor('black'), /T3D
  if (xr[0] lt 0) && (xr[1] gt 0) then plots, [0,0], yr, [zr[0],zr[0]], linestyle=2, color=truecolor('black'), /T3D
endif 
if (yr[0] lt 0) && (yr[1] gt 0) && (zr[0] lt 0) && (zr[1] gt 0) then plots, xr, [0,0], [0,0], thick=lth, linestyle=2, color=truecolor('black'), /T3D
if (xr[0] lt 0) && (xr[1] gt 0) && (zr[0] lt 0) && (zr[1] gt 0) then plots, [0,0], yr, [0,0], thick=lth, linestyle=2, color=truecolor('black'), /T3D
if (xr[0] lt 0) && (xr[1] gt 0) && (yr[0] lt 0) && (yr[1] gt 0) then plots, [0,0], [0,0], zr, thick=lth, linestyle=2, color=truecolor('black'), /T3D

; figure title
  ; get the standard character width and height
  ch_w  = !d.x_ch_size * !p.charsize  ; in device coordinates
  ch_h  = !d.y_ch_size * !p.charsize
  tmp   = convert_coord(ch_w,ch_h, /device, /to_normal)
  ch_w  = tmp[0]                      ; in normalised coordinates
  ch_h  = tmp[1]
  ; resize the title such that it fits the screen
  n_ch  = strlen(ttl)
  ttl_w = n_ch*ch_w
  if ttl_w ge 0.95 then f=0.95/ttl_w else f=1
  ttl_w = f*ttl_w
  ttl_h = f* ch_h
  xyouts, 0.5*(1.0-ttl_w),1.0-2.*ttl_h, /normal, ttl, charsize=!p.charsize*f, color=truecolor('black')

!p.charsize=oldcharsize
end
