pro nicecontour, zin, xin,yin,$
                 xrange=xrange,yrange=yrange,zrange=zrange,$
                 xtitle=xtitle, ytitle=ytitle, ztitle=ztitle,title=title,$
                 colorbar=colorbar, nlevels=nlevels, _extra=ex

; make a nice, colorful contour plot, with a colorbar (if requested)
; pseudocol is assumed to be on!
;  v1.0 mdebock, 19/08/2008

erase

x=xin
y=yin
z=zin

; get rid of NAN's
idx = where(finite(z[*,0]),cnt)
if cnt ne 0 then begin
  z = z[idx,*]
  x = x[idx]
endif
; if there are still NaNs in the other channels: replace them by zeros
idx = where(~finite(z),cnt)
if cnt ne 0 then begin
  z[idx] = 0.
endif


if ~keyword_set(xrange) then xrange=[min(x),max(x)]
if ~keyword_set(yrange) then yrange=[min(y),max(y)]
if ~keyword_set(zrange) then zrange=[min(z),max(z)]
if ~keyword_set(nlevels)then nlevels=40
if ~keyword_set(xtitle) then xtitle='X'
if ~keyword_set(ytitle) then ytitle='Y'
if ~keyword_set(ztitle) then ztitle='Z'
if ~keyword_set(title)  then title=''


; color index (from yellow to blue)
coloridx = indgen(nlevels)*245/(nlevels-1)+10

; crop the z-array to xrange and yrange
idx = where((x ge xrange[0]) AND (x le xrange[1]),cnt)
if cnt ne 0 then begin
  z=z[idx,*]
  x=x[idx]
endif
idx = where((y ge yrange[0]) AND (y le yrange[1]),cnt)
if cnt eq 0 then begin
  z=z[*,idx]
  y=y[idx]
endif

; limit the z-array to zrange
idx = where(z lt zrange[0],cnt)
if cnt ne 0 then z[idx]=zrange[0]
idx = where(z gt zrange[1],cnt)
if cnt ne 0 then z[idx]=zrange[0]



;first make the colorbar if requested
if keyword_set(colorbar) then begin
  !p.multi=2
  xbar = [0.,1.]
  ybar = [min(z) + findgen(nlevels)/(nlevels-1)*(max(z)-min(z))]
  zbar = rebin(transpose(coloridx),2,nlevels)
  contour, zbar,xbar, ybar,nlevels=nlevels,c_color=coloridx,/fill,/closed,$
           xs=1,ys=1, charsize=1.0,yticks=1,xticks=1,$
           xtickname=[' ',' '],ytickname=[' ',' '],$
           position=[0.8,0.1,0.9,0.95]
  axis,yaxis=1,ys=1,charsize=1.0,$
       ytitle=ztitle
endif
; make the contour plot
if keyword_set(colorbar) then begin
  contour, z,x,y,/fill,nlevels=nlevels,c_color=coloridx,$
           xs=1,xr=xrange,ys=1,yr=yrange,$
           charsize=1.0,$
           xtitle=xtitle, ytitle=ytitle,title=title,$
           position=[0.1,0.1,0.75,0.95], _extra=ex
  contour, z,x,y,nlevels=nlevels,$
           xs=1,xr=xrange,ys=1,yr=yrange,$
           charsize=1.0,$
           xtitle=xtitle, ytitle=ytitle,title=title,$
           position=[0.1,0.1,0.75,0.95], /noerase, _extra=ex
endif else begin
  !p.multi=0
  contour, z,x,y,/fill,nlevels=nlevels,c_color=coloridx,$
           xs=1,xr=xrange,ys=1,yr=yrange,$
           charsize=1.0,$
           xtitle=xtitle, ytitle=ytitle,title=title, _extra=ex
  contour, z,x,y,nlevels=nlevels,$
           xs=1,xr=xrange,ys=1,yr=yrange,$
           charsize=1.0,$
           xtitle=xtitle, ytitle=ytitle,title=title,$
           /noerase, _extra=ex
endelse

!p.multi=0

end