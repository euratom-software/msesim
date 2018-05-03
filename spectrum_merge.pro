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

pro spectrum_merge, inputfiles, outputfile

; routine reads in the inputfiles and merges the spectra that corresponding
; channels. The resulting spectrum is saved in the outputfile.
; The merging assumes that the channels correspond in all merged files. Therefore
; the channel array should be exactly the same in all files. Also merging
; filtered with unfiltered data is not allowed.
;
; v1.0, mdebock 23/07/2007
;
; v2.0, mdebock 15/10/2011:  * Update to use Stokes vectors as input
; v2.1, pgeelen 09/02/2012:  * Update to use 4D Stokes vectors as input 
;       mdebock                
;

;******************************************************************
;* READ THE INPUTFILES                                            *
;******************************************************************

; get the number of inputfiles
nfiles = n_elements(inputfiles)

;----------------------------------------------------------------------------------
; print some info on what's happening to the screen
;----------------------------------------------------------------------------------
print, FORMAT='("* Merging spectra from files: ")'

; loop through the inputfiles
for i=0,nfiles-1 do begin
  print, FORMAT='("        ",A)', inputfiles[i]
  ; get the inputfile
  restore, inputfiles[i]

  ; store the data in structures
  tagname = string(format='("file",(I2.2))',i)
  if i eq 0 then begin
    ; the channel arrays should be exactly the same for all files
    chanID_test    = chanID
    channels_test  = channels
    filter_test    = filterflag     
    ; central points data
    xyz0struct     = create_struct(tagname,xyz0)
    B_v0struct     = create_struct(tagname,B_v0)
    B_xyzstruct    = create_struct(tagname,B_xyz)
    B_wstruct      = create_struct(tagname,B_w)
    B_vecstruct    = create_struct(tagname,B_vec)
    C_xyzstruct    = create_struct(tagname,C_xyz)
    C_kstruct      = create_struct(tagname,C_k)
    vec0struct     = create_struct(tagname,vec0)
    Bfld0struct    = create_struct(tagname,Bfld0)
    Efld0struct    = create_struct(tagname,Efld0)
    Dshift0struct  = create_struct(tagname,Dshift0)
    Sshift0struct  = create_struct(tagname,Sshift0)
    alpha0struct   = create_struct(tagname,alpha0)
    psi0struct     = create_struct(tagname,psi0)
    lambda0struct  = create_struct(tagname,lambda0)

    ; grid point data
    gp_xyzstruct   = create_struct(tagname,gp_xyz)
    gp_velstruct   = create_struct(tagname,gp_vel)
    gp_vecstruct   = create_struct(tagname,gp_vec)
    gp_Bfldstruct  = create_struct(tagname,gp_Bfld)
    gp_Efldstruct  = create_struct(tagname,gp_Efld)
    gp_alphastruct = create_struct(tagname,gp_alpha)
    gp_psistruct   = create_struct(tagname,gp_psi)
    gp_emisstruct  = create_struct(tagname,gp_emis)

    ; resolution and emission data
    Rstruct        = create_struct(tagname,R)
    Zstruct        = create_struct(tagname,Z)
    psistruct      = create_struct(tagname,psi)
    RZpsistruct    = create_struct(tagname,RZpsi)
    Rmstruct       = create_struct(tagname,Rm)
    RZ_emisstruct  = create_struct(tagname,RZ_emis)
    R_emisstruct   = create_struct(tagname,R_emis)
    psi_emisstruct = create_struct(tagname,psi_emis)
    psi_resstruct  = create_struct(tagname,psi_res)
    R_resstruct    = create_struct(tagname,R_res)

    ; spectral data
    lambdastruct    = create_struct(tagname,lambda)
    pstokesstruct   = create_struct(tagname,pstokes)
    sstokesstruct   = create_struct(tagname,sstokes)
    stokesstruct    = create_struct(tagname,stokes)
    cwlstokesstruct = create_struct(tagname,cwlstokes)
    pibstokesstruct = create_struct(tagname,pibstokes)
    sigstokesstruct = create_struct(tagname,sigstokes)
    pirstokesstruct = create_struct(tagname,pirstokes)
  endif else begin

    ; the channel arrays should be exactly the same for all files 
    if ~array_equal(chanID_test, chanID) then begin
      print, 'ERROR: The channels of the inputfiles must be the same!'
      return
    endif
    ; and you cannot merge unfiltered with filtered data
    if  filter_test ne filterflag then begin
      print, 'ERROR: Unfiltered data cannot be merged with filtered data!'
      return
    endif

    ; central points data
    C_xyzstruct    = create_struct(C_xyzstruct, tagname,C_xyz)
    C_kstruct      = create_struct(C_kstruct, tagname,C_k)
    vec0struct     = create_struct(vec0struct, tagname,vec0)
    Dshift0struct  = create_struct(Dshift0struct, tagname,Dshift0)
    Sshift0struct  = create_struct(Sshift0struct, tagname,Sshift0)
    alpha0struct   = create_struct(alpha0struct, tagname,alpha0)

    ; grid point data
    gp_xyzstruct   = create_struct(gp_xyzstruct, tagname,gp_xyz)
    gp_velstruct   = create_struct(gp_velstruct, tagname,gp_vel)
    gp_vecstruct   = create_struct(gp_vecstruct, tagname,gp_vec)
    gp_Bfldstruct  = create_struct(gp_Bfldstruct, tagname,gp_Bfld)
    gp_Efldstruct  = create_struct(gp_Efldstruct, tagname,gp_Efld)
    gp_alphastruct = create_struct(gp_alphastruct, tagname,gp_alpha)
    gp_psistruct   = create_struct(gp_psistruct, tagname,gp_psi)
    gp_emisstruct  = create_struct(gp_emisstruct, tagname,gp_emis)

    ; resolution and emission data
    Rstruct        = create_struct(Rstruct, tagname,R)
    Zstruct        = create_struct(Zstruct, tagname,Z)
    psistruct      = create_struct(psistruct, tagname,psi)
    RZpsistruct    = create_struct(RZpsistruct, tagname,RZpsi)
    Rmstruct       = create_struct(Rmstruct, tagname,Rm)
    RZ_emisstruct  = create_struct(RZ_emisstruct, tagname,RZ_emis)
    R_emisstruct   = create_struct(R_emisstruct, tagname,R_emis)
    psi_emisstruct = create_struct(psi_emisstruct, tagname,psi_emis)
    psi_resstruct  = create_struct(psi_resstruct, tagname,psi_res)
    R_resstruct    = create_struct(R_resstruct, tagname,R_res)

    ; spectral data
    lambdastruct    = create_struct(lambdastruct, tagname,lambda)
    pstokesstruct   = create_struct(pstokesstruct, tagname,pstokes)
    sstokesstruct   = create_struct(sstokesstruct, tagname,sstokes)
    stokesstruct    = create_struct(stokesstruct, tagname,stokes)
    cwlstokesstruct = create_struct(cwlstokesstruct, tagname,cwlstokes)
    pibstokesstruct = create_struct(pibstokesstruct, tagname,pibstokes)
    sigstokesstruct = create_struct(sigstokesstruct, tagname,sigstokes)
    pirstokesstruct = create_struct(pirstokesstruct, tagname,pirstokes)
  endelse
endfor

chanID     = chanID_test
channels   = channels_test
nchan      = n_elements(channels)
filterflag = filter_test

;******************************************************************
;* MERGE CENTRAL POINT DATA                                       *
;******************************************************************
;----------------------------------------------------------------------------------
; print some info on what's happening to the screen
;----------------------------------------------------------------------------------
print, FORMAT='($,"    - Merging central point data ..................... ")'

;----------------------------------------------------------------------------------
; Get the simple merges
;----------------------------------------------------------------------------------
; for the 'zero order' information we only save the data of the first of the input files
xyz0    = xyz0struct.(0)
B_v0    = B_v0struct.(0)
B_xyz   = B_xyzstruct.(0)
B_w     = B_wstruct.(0)
B_vec   = B_vecstruct.(0)
Bfld0   = Bfld0struct.(0)
Efld0   = Efld0struct.(0)
psi0    = psi0struct.(0)
lambda0 = lambda0struct.(0)

;----------------------------------------------------------------------------------
; Get the more complex merges that required some averaging
;----------------------------------------------------------------------------------
C_xyz   = fltarr(3)
C_k     = fltarr(3)
vec0    = fltarr(3,nchan)
Dshift0 = fltarr(nchan)
Sshift0 = fltarr(nchan)
alpha0  = fltarr(nchan)
for i=0,nfiles-1 do begin
  ; get average
  C_xyz  += C_xyzstruct.(i)
  C_k    += C_kstruct.(i)
  vec0   += vec0struct.(i)
  alpha0 += alpha0struct.(i)	; we do a simple average of alpha0, because we assume 
				; that C_xyz isn't too different for the different inputfiles

  ; get maximum Doppler and Stark shift
  for j=0,nchan-1 do if abs(Dshift0struct.(i)[j]) gt abs(Dshift0[j]) then Dshift0[j] = Dshift0struct.(i)[j]
  for j=0,nchan-1 do if abs(Sshift0struct.(i)[j]) gt abs(Sshift0[j]) then Sshift0[j] = Sshift0struct.(i)[j]
endfor
C_xyz  = C_xyz/float(nfiles)
C_k    = C_k/float(nfiles)
vec0   = vec0/float(nfiles)
alpha0 = alpha0/float(nfiles)

; replace the central point structures by '0' to free up memory
xyz0struct    = 0
B_v0struct    = 0
B_xyzstruct   = 0
B_wstruct     = 0
B_vecstruct   = 0
Bfld0struct   = 0
Efld0struct   = 0
psi0struct    = 0
lambda0struct = 0
C_xyzstruct   = 0
C_kstruct     = 0
vec0struct    = 0
alpha0struct  = 0

;----------------------------------------------------------------------------------
; print some info on what's happening to the screen
;----------------------------------------------------------------------------------
print, 'done!'

;**********************************************************************************
;* MERGE GRID POINT DATA                                                          *
;**********************************************************************************

;----------------------------------------------------------------------------------
; print some info on what's happening to the screen
;----------------------------------------------------------------------------------
print, FORMAT='($,"    - Merging grid point data ........................ ")'

; the first files: 
gp_xyz   = gp_xyzstruct.(0)
gp_vel   = gp_velstruct.(0)
gp_vec   = gp_vecstruct.(0)
gp_Bfld  = gp_Bfldstruct.(0)
gp_Efld  = gp_Efldstruct.(0)
gp_alpha = gp_alphastruct.(0)
gp_psi   = gp_psistruct.(0)
gp_emis  = gp_emisstruct.(0)

; then loop through the rest of the files
for i=1,nfiles-1 do begin
  gp_xyz   = [[gp_xyz],[gp_xyzstruct.(i)]]
  gp_vel   = [[gp_vel],[gp_velstruct.(i)]]
  gp_vec   = [[gp_vec],[gp_vecstruct.(i)]]
  gp_Bfld  = [[gp_Bfld],[gp_Bfldstruct.(i)]]
  gp_Efld  = [[gp_Efld],[gp_Efldstruct.(i)]]

  gp_alpha = [gp_alpha,gp_alphastruct.(i)]
  gp_psi   = [gp_psi,gp_psistruct.(i)]
  gp_emis  = [gp_emis,gp_emisstruct.(i)]
endfor

; replace the grid point structures by '0' to free up memory
gp_xyzstruct   = 0
gp_velstruct   = 0
gp_vecstruct   = 0
gp_Bfldstruct  = 0
gp_Efldstruct  = 0
gp_alphastruct = 0
gp_psistruct   = 0
gp_emisstruct  = 0

;----------------------------------------------------------------------------------
; print some info on what's happening to the screen
;----------------------------------------------------------------------------------
print, 'done!'

;**********************************************************************************
;* MERGE RESOLUTION AND EMISSION DATA                                             *
;**********************************************************************************
;----------------------------------------------------------------------------------
; print some info on what's happening to the screen
;----------------------------------------------------------------------------------
print, FORMAT='($,"    - Merging emission and spatial resolution data ... ")'

;----------------------------------------------------------------------------------
; get minimum and maximum R and Z values and get the minimal steps in R,Z and psi
;----------------------------------------------------------------------------------
Rmax = fltarr(nfiles)
Rmin = fltarr(nfiles)
dR   = fltarr(nfiles)
Zmax = fltarr(nfiles)
Zmin = fltarr(nfiles)
dZ   = fltarr(nfiles)
dpsi = fltarr(nfiles)

for i=0,nfiles-1 do begin
  Rmax[i] = max(Rstruct.(i))
  Rmin[i] = min(Rstruct.(i))
  dR[i]   = abs(Rstruct.(i)[1]-Rstruct.(i)[0])
  Zmax[i] = max(Zstruct.(i))
  Zmin[i] = min(Zstruct.(i))
  dZ[i]   = abs(Zstruct.(i)[1]-Zstruct.(i)[0])
  dpsi[i] = abs(psistruct.(i)[1]-psistruct.(i)[0])
endfor
Rmax = max(Rmax)
Rmin = min(Rmin)
dR   = min(dR)
Zmax = max(Zmax)
Zmin = min(Zmin)
dZ   = min(dZ)
dpsi = min(dpsi)

;----------------------------------------------------------------------------------
; get overall R, Z and psi bases
;----------------------------------------------------------------------------------
nR = round((Rmax-Rmin)/dR) + 1
dR = (Rmax-Rmin)/float(nR-1)
R  = Rmin + findgen(nR)*dR

nZ = round((Zmax-Zmin)/dZ) + 1
dZ = (Zmax-Zmin)/float(nZ-1)
Z  = Zmin + findgen(nZ)*dZ

npsi = round(1.0/dpsi) + 1
dpsi = 1/float(npsi-1)
psi  = findgen(npsi)*dpsi

;----------------------------------------------------------------------------------
; get psi as a function of the new R and Z
;----------------------------------------------------------------------------------
; interpolate RZpsi to the new R and Z basis for each input file and take the average
RZpsi      = fltarr(nR,nZ)
avg_factor = fltarr(nR,nZ)
for i=0,nfiles-1 do begin
  ; the interpolation indices
  Rmaxtmp = max(Rstruct.(i))
  Rmintmp = min(Rstruct.(i))
  dRtmp   = Rstruct.(i)[1] - Rstruct.(i)[0]
  nRtmp   = n_elements(Rstruct.(i))
  Ridx   = (R - Rmintmp)/(Rmaxtmp-Rmintmp) * (nRtmp-1)

  Zmaxtmp = max(Zstruct.(i))
  Zmintmp = min(Zstruct.(i))
  dZtmp   = Zstruct.(i)[1] - Zstruct.(i)[0]
  nZtmp   = n_elements(Zstruct.(i))
  Zidx   = (Z - Zmintmp)/(Zmaxtmp-Zmintmp) * (nZtmp-1)

   ; do the interpolation of the RZpsi
  RZpsitmp = interpolate(RZpsistruct.(i), Ridx,Zidx, /grid, missing=0.0)
  idx    = where(RZpsitmp gt 0, cnt)
  if cnt ne 0 then avg_factor[idx] += 1
  ; and add up
  RZpsi += RZpsitmp
endfor
RZpsi = RZpsi/(avg_factor)


;----------------------------------------------------------------------------------
; The magnetic axis should be the same for all
;----------------------------------------------------------------------------------
Rm = Rmstruct.(0) 

;----------------------------------------------------------------------------------
; The emission intensities of each input file are interpolated to the new R, Z and 
; psi bases, making sure that the in total integrated intensity remains the same. Then
; the emission intensities of the different input files are add up
;----------------------------------------------------------------------------------
RZ_emis = fltarr(nR,nZ,nchan)
R_emis = fltarr(nR,nchan)
psi_emis = fltarr(npsi,nchan)

for i=0,nfiles-1 do begin
  ; the interpolation indices
  Rmaxtmp = max(Rstruct.(i))
  Rmintmp = min(Rstruct.(i))
  dRtmp   = Rstruct.(i)[1] - Rstruct.(i)[0]
  nRtmp   = n_elements(Rstruct.(i))
  Ridx    = (R - Rmintmp)/(Rmaxtmp-Rmintmp) * (nRtmp-1)

  Zmaxtmp = max(Zstruct.(i))
  Zmintmp = min(Zstruct.(i))
  dZtmp   = Zstruct.(i)[1] - Zstruct.(i)[0]
  nZtmp   = n_elements(Zstruct.(i))
  Zidx   = (Z - Zmintmp)/(Zmaxtmp-Zmintmp) * (nZtmp-1)

  psimaxtmp = max(psistruct.(i))
  psimintmp = min(psistruct.(i))
  dpsitmp   = psistruct.(i)[1] - psistruct.(i)[0]
  npsitmp   = n_elements(psistruct.(i))
  psiidx = (psi - psimintmp)/(psimaxtmp-psimintmp) * (npsitmp-1)

  chanidx    = indgen(nchan)

  ; do the interpolation of the emission intensity as a function of (R,Z)
  totold = total(RZ_emisstruct.(i)) * dRtmp * dZtmp
  if nchan eq 1 then begin
    RZ_emistmp = interpolate(RZ_emisstruct.(i), Ridx, Zidx, /grid, missing=0.0)
  endif else begin
    RZ_emistmp = interpolate(RZ_emisstruct.(i), Ridx, Zidx, chanidx, /grid, missing=0.0)
  endelse
  totnew = total(RZ_emistmp) * dR * dZ
  RZ_emistmp = RZ_emistmp * totold/totnew
  ; and add up
  RZ_emis += RZ_emistmp

  ; do the interpolation of the emission intensity as a function of R
  totold = total(R_emisstruct.(i)) * dRtmp
  if nchan eq 1 then begin
    R_emistmp = interpolate(R_emisstruct.(i), Ridx, /grid, missing=0.0)
  endif else begin
    R_emistmp = interpolate(R_emisstruct.(i), Ridx, chanidx, /grid, missing=0.0)
  endelse
  totnew = total(R_emistmp) * dR
  R_emistmp = R_emistmp * totold/totnew
  ; and add up
  R_emis += R_emistmp

  ; do the interpolation of the emission intensity as a function of psi
  totold = total(psi_emisstruct.(i)) * dpsitmp
  if nchan eq 1 then begin
    psi_emistmp = interpolate(psi_emisstruct.(i), psiidx, /grid, missing=0.0)
  endif else begin
    psi_emistmp = interpolate(psi_emisstruct.(i), psiidx, chanidx, /grid, missing=0.0)
  endelse
  totnew = total(psi_emistmp) * dpsi
  psi_emistmp = psi_emistmp * totold/totnew
  ; and add up
  psi_emis += psi_emistmp
endfor


;----------------------------------------------------------------------------------
; Get the new spatial resolution
;----------------------------------------------------------------------------------
resfrac = psi_resstruct.(0)[6,0]
psi_res = fltarr(7,nchan)
R_res   = fltarr(7,nchan)
for k =0,nchan-1 do begin
  ; spatial resolution as function of psi
  mx           = max(psi_emis[*,k],maxidx)
  if mx gt 0.01 then begin
    cumemis      = total(psi_emis[*,k]/mx,/cum)
    COMidx       = value_locate(cumemis, 0.5*max(cumemis))
    loweridx     = value_locate(cumemis, 0.5*(1.0-resfrac)*max(cumemis))
    upperidx     = value_locate(cumemis, 0.5*(1.0+resfrac)*max(cumemis))
    if COMidx   lt 0    then COMidx=0
    if COMidx   ge npsi then COMidx=npsi-1
    if loweridx lt 0    then loweridx=0
    if loweridx ge npsi then loweridx=npsi-1
    if upperidx lt 0    then upperidx=0
    if upperidx ge npsi then upperidx=npsi-1
    psi_res[0,k] = psi0[k]
    psi_res[1,k] = psi[maxidx]
    psi_res[2,k] = psi[COMidx]
    psi_res[3,k] = psi[loweridx]
    psi_res[4,k] = psi[upperidx]
    psi_res[5,k] = psi[upperidx]-psi[loweridx]
    psi_res[6,k] = resfrac
  endif else begin
    psi_res[0,k] = psi0[k]
    psi_res[1,k] = psi0[k]
    psi_res[2,k] = psi0[k]
    psi_res[3,k] = psi0[k]
    psi_res[4,k] = psi0[k]
    psi_res[5,k] = 0.0
    psi_res[6,k] = resfrac
  endelse

  ; spatial resolution as function of R
  mx           = max(R_emis[*,k],maxidx)
  if mx gt 0.01 then begin
    cumemis      = total(R_emis[*,k]/mx,/cum)
    COMidx       = value_locate(cumemis, 0.5*max(cumemis))
    loweridx     = value_locate(cumemis, 0.5*(1.0-resfrac)*max(cumemis))
    upperidx     = value_locate(cumemis, 0.5*(1.0+resfrac)*max(cumemis))
    if COMidx   lt 0  then COMidx=0
    if COMidx   ge nR then COMidx=nR-1
    if loweridx lt 0  then loweridx=0
    if loweridx ge nR then loweridx=nR-1
    if upperidx lt 0  then upperidx=0
    if upperidx ge nR then upperidx=nR-1
    R_res[0,k]   = norm(xyz0[0:1,k])
    R_res[1,k]   = R[maxidx]
    R_res[2,k]   = R[COMidx]
    R_res[3,k]   = R[loweridx]
    R_res[4,k]   = R[upperidx]
    R_res[5,k]   = R[upperidx]-R[loweridx]
    R_res[6,k]   = resfrac
  endif else begin
    Rtmp         = norm(xyz0[0:1,k])
    R_res[0,k]   = Rtmp
    R_res[1,k]   = Rtmp
    R_res[2,k]   = Rtmp
    R_res[3,k]   = Rtmp
    R_res[4,k]   = Rtmp
    R_res[5,k]   = 0.0
    R_res[6,k]   = resfrac
  endelse
endfor

; replace the emission structures by '0' to free up memory
Rstruct        = 0
Zstruct        = 0
psistruct      = 0
RZpsistruct    = 0
Rmstruct       = 0
RZ_emisstruct  = 0
R_emisstruct   = 0
psi_emisstruct = 0
R_resstruct    = 0
psi_resstruct  = 0

;----------------------------------------------------------------------------------
; print some info on what's happening to the screen
;----------------------------------------------------------------------------------
print, 'done!'

;**********************************************************************************
;* MERGE THE SPECTRA                                                              *
;**********************************************************************************
;----------------------------------------------------------------------------------
; print some info on what's happening to the screen
;----------------------------------------------------------------------------------
print, FORMAT='($,"    - Merging spectral data .......................... ")'

;----------------------------------------------------------------------------------
; get an overall wavelength basis
;----------------------------------------------------------------------------------
lambdamax = fltarr(nfiles)
lambdamin = fltarr(nfiles)
dlambda   = fltarr(nfiles)
for i=0,nfiles-1 do begin
  lambdamax[i] = max(lambdastruct.(i))
  lambdamin[i] = min(lambdastruct.(i))
  dlambda[i]   = lambdastruct.(i)[1] - lambdastruct.(i)[0]
endfor
lambdamax = max(lambdamax)
lambdamin = min(lambdamin)
dlambda   = min(dlambda)
nlambda   = round((lambdamax-lambdamin)/dlambda) + 1
lambda    = lambdamin + findgen(nlambda)*dlambda

;----------------------------------------------------------------------------------
; loop through the inputfiles and interpolate the Stokes vectors to 
; the overall wavelength basis (keeping the integrated intensity the same).
; Then add them up.
;----------------------------------------------------------------------------------
nstokes   = n_elements(stokesstruct.(0)[0,*,0])
pstokes   = fltarr(nlambda,nstokes,nchan)	; pi Stokes vectors (for all channels)
sstokes   = fltarr(nlambda,nstokes,nchan)	; sigma Stokes vectors (for all channels)
stokes    = fltarr(nlambda,nstokes,nchan)	; total Stokes vectors (for all channels)
for i=0,nfiles-1 do begin
  ; the interpolation indices
  lmaxtmp    = max(lambdastruct.(i))
  lmintmp    = min(lambdastruct.(i))
  dlambdatmp = lambdastruct.(0)[1] - lambdastruct.(0)[0]
  nltmp      = n_elements(lambdastruct.(i))
  lambdaidx  = (lambda - lmintmp)/(lmaxtmp-lmintmp) * (nltmp-1)
  stokesidx  = indgen(nstokes)
  chanidx    = indgen(nchan)

  ; do the interpolation of the pi-Stokes vectors
  totold = total(pstokesstruct.(i)[*,0,0]) * dlambdatmp
  if nchan eq 1 then begin
    pstokestmp = interpolate(pstokesstruct.(i), lambdaidx, stokesidx, /grid, missing=0.0)
  endif else begin
    pstokestmp = interpolate(pstokesstruct.(i), lambdaidx, stokesidx, chanidx, /grid, missing=0.0)
  endelse
  totnew = total(pstokestmp[*,0,0]) * dlambda
  pstokestmp = pstokestmp * totold/totnew
  ; and add up
  pstokes += pstokestmp

  ; do the interpolation of the sigma-Stokes vectors
  totold = total(sstokesstruct.(i)[*,0,0]) * dlambdatmp
  if nchan eq 1 then begin
    sstokestmp = interpolate(sstokesstruct.(i), lambdaidx, stokesidx, /grid, missing=0.0)
  endif else begin
    sstokestmp = interpolate(sstokesstruct.(i), lambdaidx, stokesidx, chanidx, /grid, missing=0.0)
  endelse
  totnew = total(sstokestmp[*,0,0]) * dlambda
  sstokestmp = sstokestmp * totold/totnew
  ; and add up
  sstokes += sstokestmp

  ; do the interpolation of the total Stokes vectors
  totold = total(stokesstruct.(i)[*,0,0]) * dlambdatmp
  if nchan eq 1 then begin
    stokestmp = interpolate(stokesstruct.(i), lambdaidx, stokesidx, /grid, missing=0.0)
  endif else begin
    stokestmp = interpolate(stokesstruct.(i), lambdaidx, stokesidx, chanidx, /grid, missing=0.0)
  endelse
  totnew = total(stokestmp[*,0,0]) * dlambda
  stokestmp = stokestmp * totold/totnew
  ; and add up
  stokes += stokestmp
endfor

;----------------------------------------------------------------------------------
; get the CWL stokes vector (CWL used is that of the first input file)
;----------------------------------------------------------------------------------
cwlstokes = cwlstokesstruct.(0)
for k=0,nchan-1 do begin
  lambdacwl = cwlstokes[nstokes,k]
  tmp       = min(abs(lambda-lambdacwl),idxcwl)
  cwlstokes[0:nstokes-1, k] = reform(stokes[idxcwl,*,k])
endfor

;----------------------------------------------------------------------------------
; Find the wavelength of max pi-blue, max. sigma, and max.pi-red (near the
; pi-blue, sigma and pi-red wavelenghts of the first input file)
; and the corresponding Stokes vectors
;----------------------------------------------------------------------------------
pibstokes = pibstokesstruct.(0) 
sigstokes = sigstokesstruct.(0)
pirstokes = pirstokesstruct.(0) 
for k=0,nchan-1 do begin
  ; get SN ratio
  nonzero     = where(stokes[*,0,k] gt 0.0,count)   ; where the total intensity isn't zero
  SN          = fltarr(nlambda)
  SN[nonzero] = sqrt( (stokes[nonzero,1,k]^2 + stokes[nonzero,2,k]^2)/stokes[nonzero,0,k] )
  ; approximate positions of pi and sigma
  mx = min(abs(lambda-pibstokes[4,k]), pbidx)       ; approximate position of the blue pi peak
  mx = min(abs(lambda-sigstokes[4,k]),sidx)         ; approximate position of the sigma peak
  mx = min(abs(lambda-pirstokes[4,k]), pridx)       ; approximate position of the red pi peak

  ; then find actual positions of the maxima in SN
  range    = 4.             ; find maximum with a +/-range number of pixels around the approx. positions
  ; SN maximum around blue shifted pi
  startidx = pbidx-range
  stopidx  = pbidx+range
  if startidx lt 0       then startidx=0
  if stopidx  ge nlambda then stopidx=nlambda-1
  mx       = max(SN[startidx:stopidx], pbidx)
  pbidx    = startidx + pbidx
  ; SN maximum around sigma
  startidx = sidx-range
  stopidx  = sidx+range
  if startidx lt 0       then startidx=0
  if stopidx  ge nlambda then stopidx=nlambda-1
  mx = max(SN[startidx:stopidx], sidx)
  sidx = startidx + sidx
  ; SN maximum around red shifted pi
  startidx = pridx-range
  stopidx  = pridx+range
  if startidx lt 0       then startidx=0
  if stopidx  ge nlambda then stopidx=nlambda-1
  mx = max(SN[startidx:stopidx], pridx)
  pridx = startidx + pridx
  ; pi-blue shift stokes vector
  pibstokes[0:3,k] = reform(stokes[pbidx,*,k])
  pibstokes[4,k]   = lambda[pbidx]    
  ; sigma stokes vector
  sigstokes[0:3,k] = reform(stokes[sidx,*,k])
  sigstokes[4,k]   = lambda[sidx]
  ; pi-red shift stokes vector
  pirstokes[0:3,k] = reform(stokes[pridx,*,k])
  pirstokes[4,k]   = lambda[pridx]
endfor

; replace the stokes structures by '0' to free up memory
lambdastruct    = 0
pstokesstruct   = 0
sstokesstruct   = 0
stokesstruct    = 0
cwlstokesstruct = 0
pibstokesstruct = 0
sigstokesstruct = 0
pirstokesstruct = 0

;----------------------------------------------------------------------------------
; print some info on what's happening to the screen
;----------------------------------------------------------------------------------
print, 'done!'


;**********************************************************************************
;* Save the result                                                                *
;**********************************************************************************
save, channels, chanID, xyz0, B_v0, B_xyz, B_w, B_vec, C_xyz, C_k,     $; central points data
      vec0, Bfld0, Efld0, Dshift0, Sshift0, alpha0, psi0, lambda0,     $
      gp_xyz, gp_vel, gp_vec, gp_Bfld, gp_Efld,                        $; grid point data
      gp_alpha, gp_psi, gp_emis,                                       $
      R, Z, psi, RZpsi, Rm, RZ_emis, R_emis, psi_emis, psi_res, R_res, $; resolution data
      lambda, pstokes, sstokes, stokes,                                $; spectral data
      cwlstokes, pirstokes, pibstokes, sigstokes,                      $
      filterflag,                                                      $; filterflag
      filename=outputfile                                               ; output file

end
