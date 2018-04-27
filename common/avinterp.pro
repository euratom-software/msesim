function avinterp, yold,xold,xnew
; function performs an 'averaging' interpolation:
;  - use: ynew = avinterp(yold,xold,xnew),
;
;         where: ynew(i) = mean( yold(where (xold > xnew(i-1) + xnew(i)  & xold < xnew(i) + xnew(i+1) ) )
;                                                   -------------------           -------------------
;                                                             2                              2
; v1.0 mdebock, 06-09-2007

nold = n_elements(xold)

nnew = n_elements(xnew)
ynew = fltarr(nnew)

if nnew eq 1 then begin		; if we interpolate to just 1 point the we take the weighted average of the 2 closest data points

  mn = min(abs(xold - xnew),idx)
  if (xnew-xold[idx-1]) lt (xold[idx+1]-xnew) then begin
    frac = (xold(idx)-xnew)/(xold(idx)-xold(idx-1))
    ynew = frac*yold(idx-1) + (1.0-frac)*yold(idx)
  endif else begin
    frac = (xnew-xold(idx))/(xold(idx+1)-xold(idx))
    ynew = frac*yold(idx+1) + (1.0-frac)*yold(idx)
  endelse

endif else begin		; else we use the formula given above
  for j=0,nnew-1 do begin

    if j eq 0 then begin
      xdown = (3 * xnew[j]   - xnew[j+1])/2;
      xup   = (xnew[j+1] + xnew[j])/2;
    endif else begin
      if j eq (nnew-1) then begin
        xdown = (xnew[j-1] + xnew[j])/2;
        xup   = (3 * xnew[j] - xnew[j-1])/2;
      endif else begin
        xdown = (xnew[j-1] + xnew[j])/2;
        xup   = (xnew[j+1] + xnew[j])/2;
      endelse
    endelse

    idx = where( (xold ge xdown) and (xold lt xup), count)
    if count ne 0 then begin
      ynew[j] = mean(yold[idx])
    endif else begin
      ynew[j] = interpol(yold,xold,xnew[j])
    endelse

  endfor
endelse

return, ynew

end
