; initialize post_script device 
function ps_init,name,fontsize=fontsize,aspect=aspect,idlfont=idlfont,$
                 times=times,ps=ps,portrait=portrait,zoom=zoom,over=over,$
                 plotsize=plotsize,bits_per_pixel=bits_per_pixel,_extra=ex
  
  IF NOT ARG_PRESENT(name) THEN BEGIN
    print,'Setting name to plot.eps'
    name = 'plot.eps'  
  ENDIF
  
 IF NOT KEYWORD_SET(bits_per_pixel) then bits=8 else bits= bits_per_pixel
   psys=!p
   ysys=!y
   xsys=!x
   zsys=!z
   csys=!c
   
                                ;store color table
   tvlct,r,g,b,/get
   ct = {r:r,g:g,b:b}
   
   if n_elements(plotsize) EQ 2 then begin
       plotsize=double(plotsize)
       aspect = plotsize[0]/plotsize[1]
       zoom = plotsize[0]/16.0
   endif
   IF NOT KEYWORD_SET(idlfont) THEN !p.font=0
   IF NOT KEYWORD_SET(fontsize) THEN fontsize=12
   IF NOT KEYWORD_SET(aspect) THEN aspect = (1+sqrt(5))/2.0
   IF NOT KEYWORD_SET(zoom) THEN zoom = 1.   

   if NOT KEYWORD_SET(over) then begin
       print,'Initialize PS device'
       set_plot,'PS'
       IF KEYWORD_SET(ps) THEN BEGIN
           device,/Color, File=name,portrait=portrait,encapsulated=0,bits=bits
           IF KEYWORD_SET(portrait) THEN BEGIN
               device,xsize=zoom*21.0,ysize=zoom*27.9
           ENDIF ELSE BEGIN
               device,xsize=zoom*27.9,ysize=zoom*21.0
           ENDELSE
           device,xoffset=0,yoffset=0
       ENDIF ELSE BEGIN
           device,/Color, File=name,/encapsulated,bits=bits
           device,xsize=zoom*16.0,ysize=zoom*16.0/aspect
       ENDELSE       
       device,/SYMBOL,FONT_INDEX=4
       IF KEYWORD_SET(times) THEN BEGIN
           device,/Times,/BOLD,FONT_INDEX=3,/ISOLATIN1
           device,/Times,/Italic,/ISOLATIN1
       ENDIF ELSE BEGIN
           device,/Helvetica,/BOLD,FONT_INDEX=3,/ISOLATIN1
           device,/Helvetica,/Italic,/ISOLATIN1
       ENDELSE
       device,FONT_SIZE=fontsize
   endif
   return,{p:psys,y:ysys,x:xsys,z:zsys,c:csys,ct:ct}
 end
