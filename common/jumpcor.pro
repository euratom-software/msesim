;+
; JUMPCOR, pola, r, polac, condition=condition
;
; Corrects for pi-sigma phase jumps in the MSE polarisation angle
;
; :Params:
;    pola      : input polarisation angle
;    r         : radial coordinate
;    polac     : corrected polarisation angle
; :Keywords:
;    condition : 2-element vector setting the condition:
;                      abs(pola) should be smaller than condition[1] , where r equals condition[0]
;                  AND abs(pola[i+1] - pola[i]) should always be smaller than condition[1]
;                default: condition=[0.9,!pi/4]
;-
pro jumpcor, pola,r,polac, condition=condition

  ; set the default condition
  if n_elements(condition) ne 2 then condition=[0.9,!pi/4]

  ; number of time points
  sz=size(pola)
  if sz[0] eq 1 then nt=1 else nt=sz[1] ; in case of a 1D array => assume it's just 1 radial profile

  ; initialise the array with the corrected angles
  polac=pola

  ; jump from pi to sigma
  jump      = !pi/2.0
  ; threshold for phase jump
  threshold = condition[1]/jump
  ; get the condition radius
  dummy=min(abs(r-condition[0]),iax)

  ; bring polarisation angle (which has a pi degeneracy) into [-pi/2,pi/2]
  y = pola mod !pi ; y in [-pi,+pi]
  y = y + !pi      ; y in [ 0,+2pi]
  y = y mod !pi    ; y in [ 0,  pi]
  idx = where(y ge  !pi/2., cnt)
  if cnt ne 0 then y[idx] = y[idx] - !pi ; bring [pi/2,pi] to [-pi/2, 0]
  idx = where(y lt -!pi/2., cnt)
  if cnt ne 0 then y[idx] = y[idx] + !pi ; bring [-pi,-pi/2] to [0,pi/0]

  if nt gt 1 then begin
    ; loop through the time points
    for iw=0,nt-1 do begin
      ; remove channel to channel phase jumps (twice, because you can have jumps from red-pi to blue-pi)
      y2 = unwrap(y[iw,*], jump, threshold)
      y2 = unwrap(y2, jump, threshold)
      ; apply the radius condition (twice, because you can have jumps from red-pi to blue-pi)
      if y2[iax] ge  condition[1] then y2=y2-!pi/2
      if y2[iax] ge  condition[1] then y2=y2-!pi/2
      if y2[iax] lt -condition[1] then y2=y2+!pi/2
      if y2[iax] lt -condition[1] then y2=y2+!pi/2
      ; save the corrected angle
      polac[iw,*] = y2
    endfor
  endif else begin
    ; remove channel to channel phase jumps (twice, because you can have jumps from red-pi to blue-pi)
    y2 = unwrap(y, jump, threshold)
    y2 = unwrap(y2, jump, threshold)
    ; apply the radius condition (twice, because you can have jumps from red-pi to blue-pi)
    if y2[iax] ge  condition[1] then y2=y2-!pi/2
    if y2[iax] ge  condition[1] then y2=y2-!pi/2
    if y2[iax] lt -condition[1] then y2=y2+!pi/2
    if y2[iax] lt -condition[1] then y2=y2+!pi/2
    ; save the corrected angle
    polac = y2
  endelse

end
