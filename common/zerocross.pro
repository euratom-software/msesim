function zerocross, A
; function returns the indices of A that are closest to a the zero crossing of A.
;   * this function obviously works best if A is reasonably smooth (not too noisy).
;   * it will also return the index of elements in A that are exactly zero, even if A doesn't cross zero
;   * if no zero crossings are found, -1 is returned
;
; v1.0 mdebock, 24-07-2008

signA = ceil((A)/abs(A))                               ; sign of A
idx0  = where((signA-shift(signA,1))[1:*] ne 0, count) ; indices of the points BEFORE the zero crossing,
                                                       ; be it form + to -, or vice versa!
if count eq 0 then return, -1                          ; no zero crossings found

print, idx0
A0    = abs([[A[idx0]],[A[idx0+1]]])                       ; find which one is closer to zero:
print, A0
idx0  = [[idx0],[idx0+1]]                              ; the point BEFORE the crossing   (idx0)
print, idx0
tmp   = min(A0,dimension=2,mnidx)                      ; or the point AFTER the crossing (idx0+1)
idx0  = idx0[mnidx]

return, idx0

end
