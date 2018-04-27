;-----------------------------------------------------------------------------------
; Routine: MONOTONIC
; Version: 1.10
; Author : R.Martin
; Date   : 10.12.03
;-----------------------------------------------------------------------------------;
; MONOTONIC
;
; Fuction returns true if the input array is monotonic. Note returns
; if input value has only one element.
;
; Calling Sequence
;
; result=MONOTONIC(data)
;
; result - Byte 
; data   - data type except (Pointer, Object, Structure)
;
; Optional flag: /dec - returns, true if data monotonically decreaseing
; Optional flag: /inc - returns, true if data monotonically increaseing
; Optional flag: /strict - returns, false if X(i)=X(i+1) 
;
;-----------------------------------------------------------------------------------
;

function monotonic, data, strict=strict, inc=inc, dec=dec

  type=size(data, /type)

  if (type eq 0) or (type eq 8) or (type eq 10) or (type eq 11) then return, 0B
  if (n_elements(data) le 1) then return, 0B

  if keyword_set(strict) then begin
    if keyword_set(inc) then return, array_equal(data lt data(1:*), 1B)
    if keyword_set(dec) then return, array_equal(data gt data(1:*), 1B)
    
    return, array_equal(data lt data(1:*), 1B) or $
            array_equal(data gt data(1:*), 1B)
  endif 
  if keyword_set(inc) then return, array_equal(data le data(1:*), 1B)
  if keyword_set(dec) then return, array_equal(data ge data(1:*), 1B)
    
  return, array_equal(data le data(1:*), 1B) or $
            array_equal(data ge data(1:*), 1B)

end

;--------------------------------------------------------------------------------------
; Modification History
;
; Version 1.10
;   - Correct name
;   - Add /strict option
;
