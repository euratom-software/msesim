;--------------------------------------------------
; recursive routine that goes through the old
; structure and writes string for the new structure
;--------------------------------------------------
function replace_tag, struct, tagname, data, recursive=recursive
;+
;  newstruct = replace_tag(struct, tagname, data)
;
;  Creates a new structure, which is a copy of oldStruct with the data in field tagName replaced by newData. 
;  This new data does not necessarily have to be of the same size or even type.
;  WARNING: if the recursive keyword is used and there are multiple tags with the same name on different levels,
;           then the data in  ALL tags with this name will be replaced!
;
; :Params:
;    struct   : in, required, type=structure, scalar
;              the input structure
;    tagname  : in, required, type=string, scalar
;              the tagname for which we want to replace the data. If tagname does not exist, struct is just copied!
;    data     : in, required, type=any
;              the new data for the tag
; :Keywords:
;    recursive: in, optional, type=byte, scalar
;              if set the function goes through the whole structure, replacing ALL tags with name 'tagname'.
;              if not set only the top level is searched and replaced.
; :Returns:
;    out      : copy of the input 'struct' with the data of 'tagname' replaced by 'data'
;-

  if n_params() lt 3 then begin
    print, 'USAGE: newstruct = replace_tag(struct, tagname, data)'
    return, -1
  endif

  ntags    = n_tags(struct)
  tagnames = tag_names(struct)
  for i=0,ntags-1 do begin
    if strcmp(tagnames[i],tagname,/fold_case) then begin
      newdata = data
    endif else begin
      if size(struct.(i),/type) eq 8 then begin
        if keyword_set(recursive) then begin
           newdata = replace_tag(struct.(i), tagname, data, /recursive)
        endif else begin
           newdata = struct.(i)
        endelse
      endif else begin
         newdata = struct.(i)
      endelse
    endelse
    if i eq 0 then begin
      newstruct = create_struct(tagnames[i],newdata)
    endif else begin
      newstruct = create_struct(newstruct, tagnames[i],newdata)
    endelse
  endfor

  return, newstruct
end