function  stringscan, input_string
;+
; out = STRINGSCAN( input_string )
;
; Function scans an input string for text input. It returns a cell
; array of strings.
;
; The strings in the input string can be seperated by either tabs, line
; feeds, carriage returns, ',' or ';'. Spaces are not considered to be
; a separator. Tabs, line feeds, carriage returns and ',' seperate 
; strings as elements in a row, whereas ';' begins a new row in the cell
; array. If the number of elements differs from row to row, then the row
; size (i.e. the size of the 2nd dimension) will be that of the largest
; row. All shorter rows will be filled with empty strings ''.
;
; The strings in the input_string can be enclosed in '[]', '()' or '{}',
; but they don't have to be.
;
; The text in the input strings does not have to be enclosed in quotes.
; If it is the quotes will be returned as well. As ',' and ';' are 
; used as separators, they will not be returned.
;
; E.g. the string "{ abc def, 'ghij'; klm n  , op}" will return following
; cell array: {'abc def' '''ghij'''; 'klm n' 'op'}.
;
; Leading and trailing spaces will be trimmed from the strings.
;
; The returned array has never more that 2 dimensions.
;
; :Input:
;   input_string : required, type=string, scalar
;                 string that defines an array with comma separated strings
; :Output:
;   output_array : type=string,  array
;                 the array with strings
;
;   16/08/2010 - v0.0 - Maarten De Bock (m.f.m.d.bock@tue.nl)
;                Original MATLAB version.
;   26/04/2011 - v1.0 - Maarten De Bock (m.f.m.d.bock@tue.nl)
;                Conversion to IDL.
;-

    ; Check input parameter
    if n_params() eq 0 || size(input_string, /type) ne 7 then begin
      print, ''
      print, 'USAGE:  output_array = stringscan(input_string)'
      print, '          output_array : type=string, array'
      print, '                        array with strings'
      print, '          input_string : required, type=string, scalar'
      print, '                        string that defines an array with strings.'
      print, ''
      return, -1
    end
    input_string = strtrim(input_string,2)
    if strcmp(input_string,'') then begin
      output_array=''
      return, output_array
    endif

    ; first find the rows
    rows      = strtrim(strsplit(input_string,';',/extract,/regex),2)
    n_rows    = n_elements(rows)
    ; loop through the rows and find the elements in each row
    separators = '(\{|\}|\[|\]|\(|\)|,|'+string(9b)+'|'+string(10b)+'|'+string(13b)+')';
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
    output_array = make_array(n_cols, n_rows,/string)
    c = 0
    for i=0,n_rows-1 do begin
      for j=0,n_cols-1 do begin
        if j ge element_c[i] then begin
          output_array[j,i]=''
        endif else begin
          value = elements[c]
          output_array[j,i] = value
          c = c+1
        endelse
      endfor
    endfor

    return, output_array    
end
