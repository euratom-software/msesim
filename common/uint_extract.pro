function uint_extract, A, long=long
; function extract the numerical values from the string(array) 'A' and returns an integer(array), which
; contains the unsigned integer of the numbers in the string(s) if there were numbers in the string(s),
; and returns -1 if no numerical values were found in the string(s). (max=32,767)
; if the keyword long is given than an unsigned long integer is returned (max=2,147,483,647)
; Out of range values are also returned as -1.

; get the dimension of the string array
dimstr = size(A,/dim)
if dimstr eq 0 then begin
  dimstr=1
  scalar=1
endif else scalar=0

; convert the string array into a string vector
nstr = n_elements(A)
A    = reform(A,nstr)

; convert string (array) into a byte array
b = byte(A)
; This byte array has 2 dimensions. The first dimension now the ASCII codes for each character of
; each string. The size of this dimension is the length of the largest string in the string array.
; For strings with shorter lengths the rest is filled up with zeros.

; now loop through the byte-version of each string, find the 'numerical' ASCII codes (in between 48 and 57),
; translate that in a string again and then convert it into an integer.
if keyword_set(long) then maxval=2147483647l else maxval=32767
c = lon64arr(nstr)  ; prepare the integer array
for i=0,nstr-1 do begin
  idx = where( (b[*,i] ge 48) and (b[*,i] le 57), cnt)
  if cnt ne 0 then begin
    c[i] = long64(string(b[idx,i]))
    if c[i] gt maxval then c[i]=-1
  endif else c[i]=-1
endfor
; reform c in the same dimensions as A
c = reform(c,dimstr)

if scalar eq 1 then c=c[0]

if keyword_set(long) then c = long(c) else c = fix(c)  ; prepare the integer array


return,c
end