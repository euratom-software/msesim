function str2complex, input_string, double=double
;+
; out = STR2COMPLEX( input_string )
;
; Converts an input string with a complex notation to a complex number.
; It only works with scalars
;
; The complex notation can either be:
;   str2complex('4.0 + i*3.5')
;   str2complex('4.0 + 3.5*j')
;   str2complex('4.0+i*3.5')
;   str2complex('4.0+j * 3.5')
;   ...
; Hence:
;    - the imaginary identifier can be either j or i
;    - the imaginary identifier should be multiplied with the imaginary part
;
; :Input:
;   input_string : required, type=string, scalar
;                 string that defines the complex number
; :Output:
;   out          : type=double/float, complex, scalar
;                 the complex number
; :Keywords:
;   double       : type= byte, scalar
;                 if set a double precision complex number is returned
; :History:
;   26/04/2011 - v1.0 - Maarten De Bock (m.f.m.d.bock@tue.nl)
;                first version.
;-

  ; look for the complex identifier:
  input_string = strtrim(input_string,2)
  ipos = strpos(input_string, 'i')
  if ipos eq -1 then ipos = strpos(input_string, 'j')
  if ipos eq -1 then begin
    if keyword_set(double) then out=double(input_string) else out=float(input_string)
  endif else begin
    if ipos eq 0 then begin
      tmp = 'complex(0,1)'+strmid(input_string, ipos+1)
    endif else begin
      if ipos eq strlen(input_string)-1 then begin
        tmp = strmid(input_string, 0, ipos)+'complex(0,1)'
      endif else begin
        tmp = strmid(input_string, 0, ipos)+'complex(0,1)'+strmid(input_string, ipos+1)
      endelse
    endelse
    if keyword_set(isdouble) then void = execute('out = dcomplex('+tmp+')') else void = execute('out = complex('+tmp+')')     
  endelse
   
  return, out

end