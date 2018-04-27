function read_atomdata, file, trans

; This procedure reads the atomic data file (e.g. stark_atomdata.dat)
; and returns the relative intensities for the chosen transition 'trans'
;
; input :  file     - filename of the file that contains the atomic data
;          trans    - the chosen transition: 1=H-alpha, 2=H-beta, 3=H-gamma
; returns: pisigma  - structure with 2 fields: .pi:    [4xn] array containing the linear, quadratic and
;                                                      cubic frequency shift and the intensity for the 'n' pi-lines
;                                              .sigma: [4xn] array containing the linear, quadratic and
;                                                      cubic frequency shift and the intensity for the 'n' sigma-lines
;
;  v1.0 mdebock, 25/06/2007
;

; first read in all the numerical values of the textfile into a vector:
data = read_txt(file)

; now go trough the data-array
transrd = 0		; initialise an integer that will be used to read in the type of transition
npi     = 0		; initialise an integer that will be used to read in the number of pi-lines
nsigma  = 0		; initialise an integer that will be used to read in the number of sigma-lines

n = n_elements(data)
i = 0
while (i lt n) do begin
  ; read the first integer (type of transition)
  transrd = round(data[i])
  i++

  ; read the second integer (no. of pi-lines)
  npi = round(data[i])
  i++

  ; read the pi-line data
  pi = data[i:(i+(4*npi-1))]
  pi = reform(pi,4,npi)
  i += 4*npi

  ; read the third integer (no. of sigma-lines)
  nsigma = round(data[i])
  i++

  ; read the pi-line data
  sigma = data[i:(i+(4*nsigma-1))]
  sigma = reform(sigma,4,nsigma)
  i += 4*nsigma

  ; test whether this is the transition we're looking for
  if (transrd eq trans) then break
endwhile


; for some strange reason the frequency shift factors
; assume c given in cm/s and the E-field in MV/cm. Because
; we like standard units (m/s and V/m) we need to convert
; the frequency shift factors:
pi[0,*]    *= 1e-6
sigma[0,*] *= 1e-6
pi[1,*]    *= 1e-14
sigma[1,*] *= 1e-14
pi[2,*]    *= 1e-22
sigma[2,*] *= 1e-22


; put the result in a structure
pisigma = {pi:pi, sigma:sigma}
; and return it
return, pisigma

end