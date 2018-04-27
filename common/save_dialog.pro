; similar to "dialog_filepick", but apart from the filename, it also returns the selected filter
; and automatically appends the extension of the selected filter


; event routine
pro save_dialog_event, event

  common save_dialog_info, result

  widget_control, event.top, get_uvalue=state, /no_copy

  case event.done of
    0: return
    1: begin
         result = [event.value, event.filter]
         widget_control, event.top, /destroy
       end
    2: begin
         result = ['', '']
         widget_control, event.top, /destroy
       end
  endcase
end


; main function
function save_dialog, dialog_parent=dialog_parent, file=file, path=path, title=title,$
                      overwrite_prompt=overwrite_prompt, fix_filter=fix_filter, filter=filter

  common save_dialog_info, result

  ; start the widget
  if ~keyword_set(title)            then title='Save as ...'
  if ~keyword_set(fix_filter)       then fix_filter=0
  if ~keyword_set(overwrite_prompt) then overwrite_prompt=0

  if ~keyword_set(dialog_parent) then begin
    base  = widget_base(title=title, /column, xoffset=280, yoffset=200)
  endif else begin
    base  = widget_base(title=title, /column, /modal, group_leader=dialog_parent)
  endelse
  filesel = cw_filesel(base, filenam=file, filter=filter, fix_filter=fix_filter, /save,$
                       path=path, warn_exist=overwrite_prompt)
  result = ['','']
  widget_control, base, /realize
  widget_control, base, set_uvalue=state, /no_copy
  xmanager, 'save_dialog', base

  ; get the result and add the filter-extension
  if strcmp(result[0],'') then return, result
  filelen  = strlen(result[0])
  noextlen = strpos(result[0],'.',/reverse_search)
  if (noextlen eq -1) || ((filelen-noextlen) gt 5) then begin  ; if there is no extension (extensions longer than 4 characters are
    result[0] += result[1]                                     ; considered not to be extensions), than add the filter extension
  endif else begin                                             ; else, replace the extension with the filter extension
    noext     = strmid(result[0],0,noextlen)
    result[0] = noext+result[1]
  endelse
  ; return the result
  return, result
end