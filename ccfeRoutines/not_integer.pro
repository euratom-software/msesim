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
; Function: NOT_INTEGER
; Date: 15.10.03
; Author: R.Martin
;--------------------------------------------------------------------------
; NOT_INTEGER
; 
; Function returns ture if a variable is not of type of an integer type.
; In this case an integer is any of the IDL types BYTE, INT, UINT, LONG
; ULONG
;
; Calling Sequence
;
; result=NOT_INTEGER(var)
;
; var      : any IDL-variable type
; result   : Returns FALSE(0) if var is any INTEGER-type, 
;            TRUE(1) otherwise
;
;--------------------------------------------------------------------------
;

function not_integer, var

  mask=[1,0,0,0,1,1,1,1,   1,1,1,1,0,0,0,0]
  return, byte(mask(size(var, /type)))

end
  
