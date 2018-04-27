function reset_graphic,sys,eps=eps,tif=tif, jpeg=jpeg,_extra=ex

if keyword_set(eps) then begin
   dum=ps_restore(sys,_extra=ex)
endif else begin
   if keyword_set(tif) then begin
       s = size(tif)
       if s[n_elements(s)-2] NE 7 then begin
           tif=dialog_pickfile(filter=['*.tif','*.*'],$
                               title='INIT_GRAPHIC: Save to file',$
                               default_extension='tif',/overwrite_prompt)
           if strcmp(tif,'') then tiff='plot.tif'
       endif
       img = tvrd(/TRUE)
       print,"Writing TIF file: " + tif 
       write_tiff,tif,reverse(img,3)
       dum=plot_restore(sys,/tif,_extra=ex)
   endif else if keyword_set(jpeg) then begin
       s = size(jpeg)
       if s[n_elements(s)-2] NE 7 then begin
           jpeg=dialog_pickfile(filter=['*.jpg','*.*'],$
                               title='INIT_GRAPHIC: Save to file',$
                               default_extension='jpg',/overwrite_prompt)
           if strcmp(tif,'') then tiff='plot.jpg'
       endif
       img = tvrd(/TRUE)
       print,"Writing JPEG file: " + jpeg 
       write_jpeg,jpeg,img,true=1
       dum=plot_restore(sys,/tif,_extra=ex)
   endif else dum=plot_restore(sys,_extra=ex)


endelse

return,dum

end
