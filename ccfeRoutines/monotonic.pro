;; Copyright 2018 Euratom/CCFE/TuE

;; Permission is hereby granted, free of charge, to any person obtaining a copy of this
;; software and associated documentation files (the "Software"), to deal in the Software
;; without restriction, including without limitation the rights to use, copy, modify,
;; merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
;; permit persons to whom the Software is furnished to do so, subject to the following
;; conditions:

;; The above copyright notice and this permission notice shall be included in all copies
;; or substantial portions of the Software.

;; THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
;; INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
;; PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
;; HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
;; CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
;; OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

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
