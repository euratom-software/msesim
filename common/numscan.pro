function  numscan, input_string
;+
; output_array = NUMSCAN( input_string )
;
; Function scans an input string for numeric input. It returns a numeric
; array (Double). It can handle MATLAB-style 2D arrays and complex numbers.
;
; The numbers in the input string can be seperated by either white space,
; ',' or ';'. Both whitespace and ',' seperate numbers as elements in a
; row, whereas ';' begins a new row in the array. If the number of
; elements differs from row to row, then the row size (i.e. the size of
; the 2nd dimension) will be that of the largest row. All shorter rows
; will be filled with NaNs.
; Strings that are not recognised as a number are also returned as NaN. 
;
; The numbers in the input_string can be enclosed in '[]', '()' or '{}',
; but they don't have to be.
;
; Warning: complex numbers need to be written without white space and
;          with the 'i*' in front of the imaginary part  
;          e.g. numscan('4.0+i*3.5')     returns 4.0+3.5*i
;               numscan('4.0+3.5*i')     returns 4.0
;               numscan('4.0 +3.5i')     returns [4.0, 3.5]
;               numscan('4.0 + i*3.5')   returns [4.0, NaN, 3.5*i]
;               numscan('4.0 + 3.5 * i') returns [4.0,NaN, 3.5, Nan]
;               numscan('4.0+i3.5')      returns 4.0
;
; The returned array has never more that 2 dimensions.
;
; :Input:
;   input_string : required, type=string, scalar
;                 string that defines an array with numbers
; :Output:
;   output_array : type=double/complex, array
;                 the array with numbers
;
; :History:
;   16/08/2010 - v0.0 - Maarten De Bock (m.f.m.d.bock@tue.nl)
;                Original MATLAB version.
;   26/04/2011 - v1.0 - Maarten De Bock (m.f.m.d.bock@tue.nl)
;                Conversion to IDL.
;-

    ; Check input parameter
    if n_params() eq 0 || size(input_string, /type) ne 7 then begin
      print, ''
      print, 'USAGE:  output_array = numscan(input_string)'
      print, '          output_array : type=double/complex, array'
      print, '                        array with numbers'
      print, '          input_string : required, type=string, scalar'
      print, '                        string that defines an array with numbers.'
      print, ''
      return, -1
    end
    input_string = strtrim(input_string,2)
    if strcmp(input_string,'') then begin
      output_array=!values.f_nan
      return, output_array
    endif
    ; first find the rows
    rows      = strtrim(strsplit(input_string,';',/extract,/regex),2)
    n_rows    = n_elements(rows)
    ; loop through the rows and find the elements in each row
    separators = '(\{|\}|\[|\]|\(|\)|,| |'+string(9b)+'|'+string(10b)+'|'+string(13b)+')';
    elements   = ''   ; dummy element
    element_c  = ''   ; dummy element counter
    for i=0,n_rows-1 do begin
        tmp   = strtrim(strsplit(rows[i],separators,/extract,/regex),2)
        tmp_c = 0
        for j=0,n_elements(tmp)-1 do begin
          if ~strcmp(strtrim(tmp[j],2),'') then begin
            elements =[elements,tmp[j]]
            tmp_c = tmp_c + 1
          endif
        endfor
        element_c = [element_c, tmp_c];
    endfor
    elements     = elements[1:*]
    element_c    = element_c[1:*]
    
    n_cols  = max(element_c)
    output_array = make_array(n_cols, n_rows,/dcomplex)
    c = 0
    for i=0,n_rows-1 do begin
      for j=0,n_cols-1 do begin
        if j ge element_c[i] then begin
          output_array[j,i]=complex(!values.f_nan, !values.f_nan)
        endif else begin
          tmp = elements[c]
          value = str2complex(tmp,/double)
          output_array[j,i] = value
          c = c+1
        endelse
      endfor
    endfor

    return, output_array    
end
