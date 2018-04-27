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
