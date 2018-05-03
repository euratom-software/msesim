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

;--------------------------------------------------------------------------
; Function: IS_REAL
; Date: 15.10.03
; Author: R.Martin
;--------------------------------------------------------------------------
; IS_REAL
; 
; Function returns ture if a variable is any non-Complex number type
;
; Calling Sequence
;
; result=IS_REAL(var)
;
; var      : any IDL-variable type
; result   : TRUE(1B) if var is any IDL-number other than COMPLEX
;            FALSE(0B) otherwise
;
;--------------------------------------------------------------------------
;

function is_real, var

  mask=[0,1,1,1,1,1,0,0, 0,0,0,0,1,1,1,1]
  return, byte(mask(size(var, /type)))

end
  
