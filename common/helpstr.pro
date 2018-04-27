; simple help routine for structures: writes out complete structure ...
;

;------------------------------------------
;- recursive routine that writes the info -
;------------------------------------------
pro recursive_output, str, name, level

; Some usefule info:                       type code
typearr = ['Undefined'                ,$ ; 0
           'Byte'                     ,$ ; 1
           'Integer'                  ,$ ; 2
           'Longword integer'         ,$ ; 3
           'Floating point'           ,$ ; 4
           'Double-precision floating',$ ; 5
           'Complex floating'         ,$ ; 6
           'String'                   ,$ ; 7
           'Structure'                ,$ ; 8
           'Double-precision complex' ,$ ; 9
           'Pointer'                  ,$ ; 10
           'Object reference'         ,$ ; 11
           'Unsigned Integer'         ,$ ; 12
           'Unsigned Longword Integer',$ ; 13
           '64-bit Integer'           ,$ ; 14
           'Unsigned 64-bit Integer'   ] ; 15
indent = 4 ; number of spaces to indent for each level

; if str is no structure then output what it it is
if size(str, /type) ne 8 then begin

  ; print the name and type
  space = string( format=string(format='("(A",I0,")")',indent*level),'')
  if strlen(name) gt 12 then name=string(format='(A12,"...")',name)
  print, format='($,A,A-15,":  ",A-25," = ")',space, name,typearr[size(str,/type)]
  ; construct a value
  sz = size(str) ; size vector
  if sz[0] eq 0 then begin
    if size(str,/type) eq 0 then value = '<Undefined>' else value=str
    if (size(str,/type) eq 7) && strlen(value) gt 50 then value=string(format='(A50,"...")',value)
  endif else begin
    value='Array['
    for i=1,sz[0] do begin
      if i eq sz[0] then value=string(format='(A,I0,"]")',value,sz[i]) else value=string(format='(A,I0,", ")',value,sz[i])
    endfor
  endelse
  ; and print it
  print, value

; if str is a structure then print the name and recursively call this function again for all tags in str
endif else begin
  ; print the name
  space = string( format=string(format='("(A",I0,")")',indent*level),'')
  print, format='(A,A)',space, name
  ntags    = n_tags(str)
  newnames = tag_names(str)
  newlevel = level+1
  for i=0,ntags-1 do begin
    ;recursive_output, str.(i), strlowcase(newnames[i]), newlevel
    recursive_output, str.(i), newnames[i], newlevel
  endfor
endelse

end

;----------------
;- MAIN ROUTINE -
;----------------
pro helpstr, str

; check whether the input is a structure
if size(str,/type) ne 8 then message, 'NOT a structure!'

; recursively loop through the structure and output the info
level=0
name = 'structure'
recursive_output, str, name, level

end