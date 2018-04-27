function read_txt, file, double=double

; This function reads an text spreadsheet file and returns a vector with all the numerical 
; values in the textfile. IDL-type comment-lines and empty lines are ignored.
; If the keyword 'double' is set, it return a vector of the type double, else it returns a
; float array.
;
; input :  file     - filename of the textfile that contains the data
;
; keyword: double   - if set the returned vector is of the type 'double'
;
; returns: data     - float are double vector with all the numerical values in the textfile
;
;  v1.0 mdebock, 19/07/2007
;

; first read the whole file into a string-array
text   = ''	; initialise string array
tmpstr = ''
; open the file
get_lun, funit
openr, funit, file
; read the file
while not(eof(funit)) do begin
    readf, funit, tmpstr
    tmpstr = strtrim(tmpstr,2)
    text = [text,tmpstr]
endwhile
; close the file:
free_lun, funit
close, funit

; first get rid of all the comment-lines and empty lines
datastr = ''
comment = ''
n = n_elements(text) -1
for i=1,n do begin
  if ~strlen(text[i]) then continue		; don't save empty lines
  if strcmp(text[i],';',1) then continue	; don't save comment-lines
  datastr=[datastr,text[i]]				; but do  save the rest
endfor


; now go trough the datastr-array and extract the numerical values from each string.
; 'space','tab', 'linefeed','carriage return', ';' and ',' are the accepted separators between numerical values
separators = '( |'+string(9B)+'|'+string(10B)+'|'+string(13B)+';|,)'
data = [0]			; a first dummy element
n = n_elements(datastr) -1
for i=1,n do begin
  tmpstr = strsplit(datastr[i],separators, /extract, /regex)
  tmpstr = strtrim(tmpstr,2)
  if keyword_set(double) then begin
    datatmp = double(tmpstr)
  endif else begin
    datatmp = float(tmpstr)
  endelse
  data = [data,datatmp]
endfor
; get rid of the first dummy element
data = data[1:*]

; return data
return, data

end