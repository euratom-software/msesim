;+
; result = UNWRAP(y, jump, threshold)
;
; Unwraps phase jumps. (G. Cunningham)
;
; :Params:
;    y         : input signlapolarisation angle
;    jump      : corrective jump to make when the next point is over the threshold
;    threshold : threshold (in fractions of the jump) that decides whether the next point did have a phase jump are not
;-
function unwrap, y, jump, threshold

  nend=n_elements(y) - 1
  yu=y
  n=where(diff(yu) gt jump * threshold)

  if n[0] ge 0 then begin
    for nn=0, n_elements(n)-1 do begin
      yu[n[nn]+1 : nend] = yu[n[nn]+1 : nend] - jump
    endfor
  endif

  n=where(diff(yu) lt -(jump * threshold))
  if n[0] ge 0 then begin
    for nn=0, n_elements(n)-1 do begin
      yu[n[nn]+1 : nend] = yu[n[nn]+1 : nend] + jump
    endfor
  endif

  return, yu

end
