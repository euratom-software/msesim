PRO plotzoom_event, ev
WIDGET_CONTROL, ev.TOP, get_uvalue=uval
plotx=uval.x
ploty=uval.y
optstr=uval.optionstring


if ev.ID eq ev.TOP then begin
	
	Widget_Control, uval.drawid, XSize=ev.x, YSize=ev.y

	WIDGET_CONTROL, uval.drawid, get_value=draw
	WSET, draw
	
	dum=execute(strjoin(['plot,[2*uval.yrange(1),2*uval.yrange(1)],xrange=uval.xrange,yrange=uval.yrange,xstyle=1,ystyle=1',optstr]))
	if uval.nplots eq 1 then oplot, *plotx, *ploty, color=uval.colors
	if uval.nplots gt 1 then begin
		for i=0,uval.nplots-1 do begin
			oplot, (*plotx)(*,i),(*ploty)(*,i),color=uval.colors(i)
		endfor
	endif




endif else begin

if ev.press eq 1 then begin
	
	x1=convert_coord(ev.x,ev.y,/device,/to_data)	
	uval.xrange(0)=x1(0)
	uval.yrange(0)=x1(1)
	
	WIDGET_CONTROL, ev.TOP, set_uvalue=uval
endif


if ev.release eq 1 then begin
	
	x2=convert_coord(ev.x,ev.y,/device,/to_data)
	uval.xrange(1)=x2(0)
	uval.yrange(1)=x2(1)

uval.xrange=uval.xrange[sort(uval.xrange)]
uval.yrange=uval.yrange[sort(uval.yrange)]
if uval.xrange(1)-uval.xrange(0) eq 0 then begin
uval.xrange=uval.xro
uval.yrange=uval.yro
endif


WIDGET_CONTROL, uval.drawid, get_value=draw
WSET, draw


	dum=execute(strjoin(['plot,[2*uval.yrange(1),2*uval.yrange(1)],xrange=uval.xrange,yrange=uval.yrange,xstyle=1,ystyle=1',optstr]))
	if uval.nplots eq 1 then oplot, *plotx, *ploty, color=uval.colors
	if uval.nplots gt 1 then begin
		for i=0,uval.nplots-1 do begin
			oplot, (*plotx)(*,i),(*ploty)(*,i),color=uval.colors(i)
		endfor
	endif

endif

if (ev.release EQ 2) then begin
	
	x=convert_coord(ev.x,ev.y,/device,/to_data)
	xspan=uval.xrange(1)-uval.xrange(0)
	yspan=uval.yrange(1)-uval.yrange(0)
	uval.xrange=[x(0)-0.75*xspan,x(0)+0.75*xspan]
	uval.yrange=[x(1)-0.75*yspan,x(1)+0.75*yspan]
	WIDGET_CONTROL, uval.drawid, get_value=draw
	WSET, draw
	

	dum=execute(strjoin(['plot,[2*uval.yrange(1),2*uval.yrange(1)],xrange=uval.xrange,yrange=uval.yrange,xstyle=1,ystyle=1',optstr]))
	if uval.nplots eq 1 then oplot, *plotx, *ploty, color=uval.colors
	if uval.nplots gt 1 then begin
		for i=0,uval.nplots-1 do begin
			oplot, (*plotx)(*,i),(*ploty)(*,i),color=uval.colors(i)
		endfor
	endif

endif







if ev.press eq 4 then begin
	uval.pantemp=convert_coord(ev.x,ev.y,/device,/to_data)	
	WIDGET_CONTROL, ev.TOP, set_uvalue=uval
endif


if ev.release eq 4 then begin
	
	pan2=convert_coord(ev.x,ev.y,/device,/to_data)
	deltax=pan2(0)-uval.pantemp(0)
	deltay=pan2(1)-uval.pantemp(1)

if (deltax eq 0) and (deltay eq 0) then begin
plotzoom, *plotx, *ploty,colors=uval.colors,optionstring=optstr
endif else begin
		uval.xrange=uval.xrange-deltax
		uval.yrange=uval.yrange-deltay


		WIDGET_CONTROL, uval.drawid, get_value=draw
		WSET, draw


		dum=execute(strjoin(['plot,[2*uval.yrange(1),2*uval.yrange(1)],xrange=uval.xrange,yrange=uval.yrange,xstyle=1,ystyle=1',optstr]))
		if uval.nplots eq 1 then oplot, *plotx, *ploty, color=uval.colors
		if uval.nplots gt 1 then begin
			for i=0,uval.nplots-1 do begin
				oplot, (*plotx)(*,i),(*ploty)(*,i),color=uval.colors(i)
			endfor
		endif

endelse
endif












;if (ev.release EQ 2) then begin
;	plotzoom, *plotx, *ploty,colors=uval.colors,optionstring=optstr
;endif



endelse

WIDGET_CONTROL, ev.TOP, set_uvalue=uval
END




PRO plotzoom_cleanup, base
WIDGET_CONTROL, base, get_uvalue=uval
PTR_FREE, uval.x
PTR_FREE, uval.y
END




PRO plotzoom,x,y,optionstring=optionstring,colors=colors

dim=size(x)
if dim(0) eq 1 then nplots=1 else nplots=dim(2)


if n_elements(y) eq 0 then begin
	y=x
	if dim(0) eq 1 then x=indgen(n_elements(x),/double)
	if dim(0) eq 2 then begin
		for i=0,nplots-1 do begin
			x(*,i)=indgen(dim(1),/double)
		endfor
	endif
endif


if (not keyword_set(optionstring)) then optionstring=''
if (not keyword_set(colors)) then begin
	colors=intarr(nplots)
	for i=0,nplots-1 do begin
	colors(i)=!P.color
	endfor
endif

base=WIDGET_BASE(/COLUMN, tlb_size_events=1,title='Plotzoom')
drawid=WIDGET_DRAW(base, xsize=400,ysize=400, /button_events)


WIDGET_CONTROL, base, /REALIZE
WIDGET_CONTROL, drawid, get_value=draw
WSET, draw
xrange=[min(x),max(x)]
yrange=[min(y),max(y)]
xro=xrange ;original xrange
yro=yrange
dum=execute(strjoin(['plot,[2*yrange(1),2*yrange(1)],xrange=xrange,yrange=yrange,xstyle=1,ystyle=1',optionstring]))
if nplots eq 1 then oplot, x, y, color=colors
if nplots gt 1 then begin
	for i=0,nplots-1 do begin
		oplot, x(*,i),y(*,i),color=colors(i)
	endfor
endif

pantemp=fltarr(3)
uval={x:ptr_new(x),y:ptr_new(y),optionstring:optionstring,nplots:nplots,colors:colors,drawid:drawid,xrange:xrange,yrange:yrange,xro:xro,yro:yro,pantemp:pantemp}
WIDGET_CONTROL, base, set_uvalue=uval
XMANAGER, 'plotzoom', base, cleanup='plotzoom_cleanup'

END
