;--------------------------------------------------------------------------
; Function: IS_INTEGER
; Date: 15.10.03
; Author: R.Martin
;--------------------------------------------------------------------------
; EXISTS
; 
; Checks if a variable exists or not.

; Calling Sequence
;
; result=EXISTS(var)
;
; var      : any IDL-variable type
; result   : Returns TRUE(1B) if var exists
;            FALSE(1) otherwise
;
;--------------------------------------------------------------------------
;

function exists, var

  return, (n_elements(var) ne 0)

end
  
