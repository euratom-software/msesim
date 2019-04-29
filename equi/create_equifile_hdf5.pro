pro create_equifile, shotno=shotno, time=time, hdfname=hdfname, filename=filename
;
; This function uses the 'standard' IDL function 'read_flux' to get
; in the B-field and normalised fluxcoordinates for the given shotno.
; It then gets the B-field, the magnetic axis and normalised fluxcoordinates 
; at the given time, reforms the B-field-array such that B[*,Ridx,Zidx]=[BR,BZ,BPhi] and
; saves it to the file 'filename'

; check if keywords are given correctly
if not keyword_set(shotno) then begin
  print, 'ERROR: You must specify a shotnumber!'
  return
endif
if not keyword_set(time) then begin
  print, 'ERROR: You must specify a time!'
  return
endif
if not keyword_set(filename) then begin
  print, 'ERROR: You must specify a filename!'
  return
endif

; Read in the data
data = read_flux(shotno,/nofc) ; only works on MAST

; get the index of the time we're interested in
mnt = min(abs(data.taxis.vector-time), tidx)

; get the major radius (m)
R  = data.xaxis.vector
nR = n_elements(R)
; get the Z axis (m)
Z = data.yaxis.vector
nZ = n_elements(Z)

; get the B-field vectors at the right time and transform them from
; B[Ridx,Zidx,tidx,*]=[BR,BPhi,BZ] into B[*,Ridx,Zidx]=[BR,BZ,BPhi]
Bfldtmp     = transpose(reform(data.bfield.vector[*,*,tidx,*]),[2,0,1])
Bfld        = fltarr(3,nR,nZ)
Bfld[0,*,*] = Bfldtmp[0,*,*]
Bfld[1,*,*] = Bfldtmp[2,*,*]
Bfld[2,*,*] = Bfldtmp[1,*,*]

; get the normalised flux coordinates at the right time
fluxcoord = data.fluxcoordinates.psin[*,*,tidx]

; get the magnetic axis at the right time
Rm        = data.magnetic_axis.position[0,tidx]

; save this into the given file
save, R,Z,Bfld,fluxcoord,Rm, filename=filename


end
