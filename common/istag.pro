function istag, strct, tag,loc=loc
; return = ISTAG(STRCT, TAG, LOC=LOC)
;
;    STRCT: structure to check
;    TAG  : string with tag name
;    LOC  : index of the tag (if found) : strct.tag == strct.(loc)
;
; function checks whether 'tag' is a tag of the structure 'strct'. It returns 1 if it is,
; and returns 0 if it isn't. The keyword 'loc' returns the index of 'tag' within 'strct'.
;

; list of tag names in 'strct'
tags=tag_names(strct)
ntags=n_elements(tags)
; loop through the tag names and check whether 'tag' is one of them
; and return
for i=0,ntags-1 do begin
  if strcmp(tags[i],tag,/fold_case) then begin
    loc=i
    return, 1
  endif
endfor

; tag not found:
return, 0

end
