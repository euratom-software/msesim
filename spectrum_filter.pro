pro spectrum_filter, settingfile, inputfile, outputfile

; routine filters the (raw) spectrum in the inputfile, based on the
; filter parameters in the settingsfile, and saves the filtered 
; spectrum in the outputfile
;
; v1.0, mdebock 20/07/2007
;
; v1.1, daussems 29/09/2011: * Windows/Unix compatibility ensured
;       mdebock              * Position optimal pi-red, sigma and pi-blue and corresponding Stokes vectors determined 


;----------------------------------------------------------------------------------
; Distinguish between Windows and UNIX directory separators
;----------------------------------------------------------------------------------
if strcmp( !version.os,'Win32',/fold_case) then sep='\' else sep='/'

;**********************************************************************************
;* READ FILTER SETTINGS AND INPUTFILE                                             *
;**********************************************************************************

;----------------------------------------------------------------------------------
; filter settings
;----------------------------------------------------------------------------------
input = read_setting(settingfile)

fltfile = input.filter.filterfile
if ~strcmp(fltfile,'none',/fold_case) then fltfile= 'filter'+sep+fltfile
mapfile = input.filter.mappingfile
if ~strcmp(mapfile,'none',/fold_case) then mapfile= 'filter'+sep+mapfile
fltshape = input.filter.filtershape
fltwidth = input.filter.filterwidth
fltangle = input.filter.angle *!dtor
fltneff  = input.filter.neff

;----------------------------------------------------------------------------------
; Read the inputdata
;----------------------------------------------------------------------------------
restore, inputfile

; get some inportant data
nchan   = n_elements(channels)			; number of channels
dlambda = lambda[1]-lambda[0]			; step in wavelength
nlambda = n_elements(lambda)			; number of wavelengths


;**********************************************************************************
;* Apply the filter                                                               *
;**********************************************************************************

;----------------------------------------------------------------------------------
; print some info on what's happening to the screen
;----------------------------------------------------------------------------------
print,FORMAT='("* Applying filter:")'
print,FORMAT='("  - Correction for finite spread in incident angle applied: half angle=",F4.2,"deg. , n_eff=",F4.2)',$
      fltangle*!radeg, fltneff

;----------------------------------------------------------------------------------
; get the filter-kernel
;----------------------------------------------------------------------------------
if (fltfile eq 'none') then begin

; REMINDER: NEED TO MAKE THE FILTER CWL COINCIDE WITH THE SPECTRUM CWL FOR EACH CHANNEL
  flttrans = fltarr(nlambda,nchan)
  fwhm     = fltarr(nchan) + fltwidth
  trans    = fltarr(nchan) + 1.0
  cwl      = reform(cwlstokes[4,*])
  wlshift  = fltarr(nchan)
  cwlnew   = fltarr(nchan)
  fwhmnew  = fltarr(nchan)
  transnew = fltarr(nchan)
  for k=0,nchan-1 do begin
    ; => same filter function for every channel/fibre bundle, hence a 1D vector
    flt = calc_filter(fltshape, fltwidth, cwl[k], lambda)
    ; apply the effect of the finite angle spread on the filter
    fltstruct     = angle_filter(flt, lambda, fltangle, fltneff)
    flttrans[*,k] = fltstruct.fltrc            ; use the recentred filter function (i.e. the
                                               ; the max. transmission is still at the CWL)
    wlshift[k]    = fltstruct.wlshift          ; the wavelength shift we would expect (but we're not using it)
    cwlnew[k]     = fltstruct.cwl              ; the new CWL we would expect (but we're not using it)
    fwhmnew[k]    = fltstruct.fwhm
    transnew[k]   = fltstruct.T
  endfor
endif else begin
  ; => different filter function for every channel/fibre bundle, hence a 2D array
  fltstruct = read_filter(fltfile, mapfile, chanID)
  flttrans  = fltstruct.flt	; filter transmission function
  fltlambda = fltstruct.lambda  ; wavlength basis of the measured filters
  cwl       = fltstruct.cwl     ; central wavelengths of the measured filters
  fwhm      = fltstruct.fwhm    ; FWHM of the measured filters
  trans     = fltstruct.trans   ; peak transmission of the measured filters

  ; resample the filter wavelength vector 'fltlambda' to the spectrum wavelength vector
  lambdaidx = (lambda-min(fltlambda))/(max(fltlambda)-min(fltlambda)) * (n_elements(fltlambda)-1)
  chanidx   = indgen(nchan)
  if nchan gt 1 then flttrans = interpolate(flttrans, lambdaidx, chanidx, /grid) else flttrans = interpolate(flttrans, lambdaidx)


  ; apply the effect of the finite angle spread on the filters
  ; (we assume the filter measurement was done with a pinhole illumination of 1mm
  ;  and a 75mm focal length lens
  ;    <=> 0.38 degrees half angle of the spread in incidence angles
  ;    <=> negligable shift (-0.05A) and no increase in width (0.00A) 
  ;  so effectively a it is a 'perfect' measurement)
  wlshift  = fltarr(nchan)
  cwlnew   = fltarr(nchan)
  fwhmnew  = fltarr(nchan)
  transnew = fltarr(nchan)
  for k=0,nchan-1 do begin
    fltstruct = angle_filter(flttrans[*,k], lambda, fltangle, fltneff)
    flttrans[*,k] = fltstruct.flt            ; use the shifted filter function
    wlshift[k]    = fltstruct.wlshift
    cwlnew[k]     = fltstruct.cwl
    fwhmnew[k]    = fltstruct.fwhm
    transnew[k]   = fltstruct.T
  endfor

endelse


;----------------------------------------------------------------------------------
; loop through the channels, and for each channel through the (nonzero-intensity)
; wavelengths
;----------------------------------------------------------------------------------
for k=0,nchan-1 do begin

  ;----------------------------------------------------------------------------------
  ; print some info on what's happening to the screen
  ;----------------------------------------------------------------------------------
  print,FORMAT='("  - Filtering spectrum ",A," at R=",F4.2,"m")', chanID[k], R_res[2,k]
  print,FORMAT='("      original filter specifications:")'
  print,FORMAT='("         CWL               = ",F7.2,"A")', cwl[k]
  print,FORMAT='("         FWHM              = ",F7.2,"A")', fwhm[k]
  print,FORMAT='("         peak Transmission = ",F7.2,"%")', trans[k]*100.
  print,FORMAT='("      after finite incident angle correction:")'
  print,FORMAT='("         CWL               = ",F7.2,"A (shift = ",F5.2,"A)")', cwlnew[k], wlshift[k]
  if (fltfile eq 'none') then begin
    print,FORMAT='("           (this new CWL is not used, because of the filter was calculated not measured!)")'
  endif
  print,FORMAT='("         FWHM              = ",F7.2,"A")', fwhmnew[k]
  print,FORMAT='("         peak Transmission = ",F7.2,"%")', transnew[k]*100.

  ;----------------------------------------------------------------------------------
  ; make some temporary stokes vectors
  ;----------------------------------------------------------------------------------
  pstokes_tmp = fltarr(nlambda,4)
  sstokes_tmp = fltarr(nlambda,4)


  ;----------------------------------------------------------------------------------
  ; Find the index of the CWL of the filter in the lambda vector
  ;----------------------------------------------------------------------------------
  cflt = value_locate(lambda,cwlnew[k])
  if cflt lt 0 then cflt=0
  if cflt gt nlambda-1 then cflt=nlambda-1

  ;----------------------------------------------------------------------------------
  ; Sweep the filter kernel over the spectrum
  ; (i.e. the CWL goes from one edge to the other edge of the spectrum)
  ;----------------------------------------------------------------------------------
  for l=0, nlambda-1 do begin
    ;----------------------------------------------------------------------------------
    ; truncate the edges if necessary
    ;----------------------------------------------------------------------------------
    if (l lt cflt) then begin
      fltmin = cflt-l
      lmin   = 0
    endif else begin
      fltmin = 0
      lmin   = l-cflt
    endelse
    if ( l gt cflt ) then begin
      fltmax = nlambda+cflt-l-1
      lmax   = nlambda-1
    endif else begin
      fltmax = nlambda-1
      lmax   = nlambda-cflt+l-1
    endelse

    ;----------------------------------------------------------------------------------
    ; calculate the value of the stokes vectors for the filter CWL at lambda[l]
    ;----------------------------------------------------------------------------------
    pstokes_tmp[l,0] = total(flttrans[fltmin:fltmax,k] * pstokes[lmin:lmax,0,k])*dlambda
    pstokes_tmp[l,1] = total(flttrans[fltmin:fltmax,k] * pstokes[lmin:lmax,1,k])*dlambda
    pstokes_tmp[l,2] = total(flttrans[fltmin:fltmax,k] * pstokes[lmin:lmax,2,k])*dlambda
    pstokes_tmp[l,3] = total(flttrans[fltmin:fltmax,k] * pstokes[lmin:lmax,3,k])*dlambda

    sstokes_tmp[l,0] = total(flttrans[fltmin:fltmax,k] * sstokes[lmin:lmax,0,k])*dlambda
    sstokes_tmp[l,1] = total(flttrans[fltmin:fltmax,k] * sstokes[lmin:lmax,1,k])*dlambda
    sstokes_tmp[l,2] = total(flttrans[fltmin:fltmax,k] * sstokes[lmin:lmax,2,k])*dlambda
    sstokes_tmp[l,3] = total(flttrans[fltmin:fltmax,k] * sstokes[lmin:lmax,3,k])*dlambda
    
    ;----------------------------------------------------------------------------------
    ; calculate the value of the stokes vectors for the filter CWL at actual CWL
    ;----------------------------------------------------------------------------------
    if l eq cflt then begin
       cwlstokes[0:3,k] = pstokes_tmp[l,*] + sstokes_tmp[l,*]
       cwlstokes[4,k]   = lambda[l]
    endif
        
  endfor

  ;----------------------------------------------------------------------------------
  ; Save the temporary stokes vectors
  ;----------------------------------------------------------------------------------
  pstokes[*,*,k] = pstokes_tmp
  sstokes[*,*,k] = sstokes_tmp
  stokes[*,*,k]  = pstokes_tmp + sstokes_tmp

  ;----------------------------------------------------------------------------------
  ; Find the wavelength of max pi-blue, max. sigma, and max.pi-red
  ; and the corresponding Stokes vectors
  ;----------------------------------------------------------------------------------    
  ; get SN ratio
  nonzero     = where(stokes[*,0,k] gt 0.0,count)   ; where the total intensity isn't zero
  SN          = fltarr(nlambda)
  SN[nonzero] = sqrt( (stokes[nonzero,1,k]^2 + stokes[nonzero,2,k]^2)/stokes[nonzero,0,k] )

  ; positions of pi and sigma before filtering
  mx = min(abs(lambda-pibstokes[4,k]), pbidx)       ; approximate position of the blue pi peak
  mx = min(abs(lambda-sigstokes[4,k]),sidx)         ; approximate position of the sigma peak
  mx = min(abs(lambda-pirstokes[4,k]), pridx)       ; approximate position of the red pi peak
  ; range of in which to find the pi and sigma position <=> 1/4 the distance between pi-blue and pi-red
  range = round(0.25*(pridx-pbidx))  
  
  ; then find actual positions of the maxima in SN
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

  print,FORMAT='(%"\b\b\b\b","done!")'
endfor


;set the filterflag to 1 (because the spectrum is now filtered)
filterflag=1

;-----------------------------------------------------------------------------------------------------
; Save the result
;-----------------------------------------------------------------------------------------------------
save, channels, chanID, xyz0, B_v0, B_xyz, B_w, B_vec, C_xyz, C_k,     $; central points data
      vec0, Bfld0, Efld0, Dshift0, Sshift0, alpha0, psi0, lambda0,     $
      gp_xyz, gp_vel, gp_vec, gp_Bfld, gp_Efld,                        $; grid point data
      gp_alpha, gp_psi, gp_emis,                                       $
      R, Z, psi, RZpsi, Rm, RZ_emis, R_emis, psi_emis, psi_res, R_res, $; resolution data
      lambda, pstokes, sstokes, stokes,                                $; spectral data
      cwlstokes, pirstokes, pibstokes, sigstokes,                      $
      filterflag,                                                      $; filterflag
      filename=outputfile; output file
end
