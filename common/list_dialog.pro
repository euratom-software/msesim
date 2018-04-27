; shows a dialog window that asks the users to pick one of the possible options out of a list
; is returns the index of the selected element or -1 if cancel was clicked

; event routine
pro list_dialog_event, event

  common list_dialog_info, ok, cancel, listid, listidx

  id = event.id

  if id eq ok then begin
    listidx = widget_info(listid,/list_select)
    widget_control, event.top, /destroy
  endif
  if id eq cancel then begin
    listidx = -1
    widget_control, event.top, /destroy
  endif
end

; main function
function list_dialog, dialog_parent=dialog_parent, text=text, list=list, default=default, title=title

  common list_dialog_info, ok, cancel, listid, listidx

  ; start the widget
  if ~keyword_set(title)         then title='Select dialog'
  if ~keyword_set(text)          then text='Select one of the following ...'
  if ~keyword_set(default)       then default=0

  if ~keyword_set(dialog_parent) then begin
    base  = widget_base(title=title, /column, xoffset=280, yoffset=200)
  endif else begin
    base  = widget_base(title=title, /column, /modal, group_leader=dialog_parent, align_center=keyword_set(center))
  endelse
  for i=0,n_elements(text)-1 do begin
    lbl     = widget_label(base, value=text[i], /align_left)
  endfor
  if n_elements(list) lt 4 then ysize=n_elements(list) else ysize=4
  listid  = widget_list(base, value=list, ysize=ysize)
  widget_control, listid, set_list_select=default
  row     = widget_base(base, /row)
  ok      = widget_button(row, value='OK')
  cancel  = widget_button(row, value='Cancel')


  widget_control, base, /realize
  xmanager, 'list_dialog', base

  return, listidx
end