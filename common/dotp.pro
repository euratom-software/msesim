function dotp, V1, V2
; function returns the dot product of 2 vectors
; e.g.: ad+bf+cg = dotp([a,b,c],[d,f,g])
  return, total(V1*V2)
end