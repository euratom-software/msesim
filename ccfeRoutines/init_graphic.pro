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

function init_graphic,eps=eps,sys=sys,_extra=ex
if not keyword_set(sys) then begin
    if keyword_set(eps) then begin
        s = size(eps)
        if s[n_elements(s)-2] NE 7 then begin
            eps=dialog_pickfile(filter=['*.eps','*.*'],$
                                title='INIT_GRAPHIC: Save to file',$
                                default_extension='eps',/overwrite_prompt)
            if strcmp(eps,'') then eps='plot.eps'
        endif   
        sys=ps_init(eps,_extra=ex)
    endif else begin
        sys=plot_init(_extra=ex)
    endelse
endif 
    
return,sys

end
