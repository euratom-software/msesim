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
