; initialize post_script device 
function plot_init,fontsize=fontsize,aspect=aspect,times=times,window=window,$
                   zoom=zoom,ps=ps,portrait=portrait,tif=tif,jpeg=jpeg,$
                   plotsize=plotsize,_EXTRA=ex
  
  
   psys=!p
   ysys=!y
   xsys=!x
   zsys=!z
   csys=!c
   x_chs=!d.x_ch_size
   y_chs=!d.y_ch_size
                                ;store color table
   tvlct,r,g,b,/get
   ct = {r:r,g:g,b:b}
   
   if n_elements(plotsize) EQ 2 then begin
       plotsize=double(plotsize)
       aspect = plotsize[0]/plotsize[1]
       zoom = plotsize[0]/16.0
   endif
   IF NOT KEYWORD_SET(fontsize) THEN fontsize=12
   IF NOT KEYWORD_SET(aspect) THEN aspect = (1+sqrt(5))/2.0
   if NOT KEYWORD_SET(zoom) THEN zoom=1.
   IF KEYWORD_SET(tif) or keyword_set(jpeg) THEN BEGIN
       pixmap=1b 
       scal = 150./(2.54*16.0) ; 150 dpi
       !p.thick = round(zoom*scal)
       !x.thick = round(zoom*scal)
       !y.thick = round(zoom*scal)
       !z.thick = round(zoom*scal)
       !p.charthick = round(zoom*scal)
   ENDIF ELSE BEGIN 
       pixmap=0b
       scal = 1
   ENDELSE

   set_plot,'X'   
   if N_ELEMENTS(window) EQ 0 then begin
       if keyword_set(ps) then begin
           if keyword_set(portrait) then begin
               window,/free,xsize=scal*zoom*21.0*!D.X_PX_CM,$
                      ysize=scal*zoom*27.9*!D.X_PX_CM,$
                      title=name,pixmap=pixmap
           endif else begin
               window,/free,xsize=scal*zoom*27.9*!D.X_PX_CM,$
                      ysize=scal*zoom*21.0*!D.X_PX_CM,$
                      title=name,pixmap=pixmap
           endelse
       endif else begin
           window,/free,xsize=scal*zoom*16.0*!D.X_PX_CM,$
                  ysize=scal*zoom*16.0*!D.X_PX_CM/aspect,$
                  title=name,pixmap=pixmap
       endelse
   endif else begin
       if !D.window EQ -1 then begin
           if (window GE 0) THEN BEGIN 
               window,window,xsize=scal*zoom*16.0*!D.X_PX_CM,$
                      ysize=scal*zoom*16.0*!D.X_PX_CM/aspect,$
                      title=name,pixmap=pixmap
           endif else if window EQ 0 then begin
               window,0,xsize=scal*zoom*16.0*!D.X_PX_CM,$
                      ysize=scal*zoom*16.0*!D.X_PX_CM/aspect,$
                      title=name,pixmap=pixmap
           endif
       endif else begin
           if (window GE 0) AND (NOT KEYWORD_SET(tif)) then wset,window
       endelse
   endelse
   if not keyword_set(idlfonts) then begin
                                ; calculate the 12 point Postscript charsize
     pt=scal*fontsize*0.0353*!D.X_PX_CM
     charsize=[!D.X_CH_SIZE,!D.Y_CH_SIZE]
     device,set_character_size=[6*pt/10,pt]
   endif
   IF KEYWORD_SET(times) THEN BEGIN
     !p.font=-6
   ENDIF ELSE BEGIN
     !p.font=-3
   ENDELSE
   
   return,{p:psys,y:ysys,x:xsys,z:zsys,c:csys,ct:ct}
 end
