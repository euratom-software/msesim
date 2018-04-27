;***************************************************************************************
;* Function reads in the SPE-file format of the PI Acton WinSpec spectrometer software *
;***************************************************************************************

; routine to get data (double, char or byte) from a binary file
;---------------------------------------------------------------
pro getd,lun,pt,tp,var,rep=rep

; lun = 'logical unit number' <=> file identifier
; pt  = pointer to the location in the file where the data is found
; tp  = data type: 'float', 'double', 'byte', 'int', 'uint', 'long', 'ulong',
;                  'char', 'short', 'word', 'dword'
;       (Last 4 data types are described in WinSpec manual:
;        char=signed byte, short=int, word=uint, dword=ulong)
; var = the variable to store the data to
; keyword rep = makes var into an array of length 'rep'

; set up the data type of the variable
if strcmp(tp,'float',/fold_case) then var=0.0
if strcmp(tp,'double',/fold_case) then var=0.0D
if strcmp(tp,'byte',/fold_case) then var=0B
if strcmp(tp,'int',/fold_case) then var=0
if strcmp(tp,'uint',/fold_case) then var=0U
if strcmp(tp,'uint',/fold_case) then var=0U
if strcmp(tp,'long',/fold_case) then var=0L
if strcmp(tp,'ulong',/fold_case) then var=0UL
if strcmp(tp,'char',/fold_case) then var=0B
if strcmp(tp,'short',/fold_case) then var=0
if strcmp(tp,'word',/fold_case) then var=0U
if strcmp(tp,'dword',/fold_case) then var=0UL

; turn the variable into an array of length 'rep' if the rep-keyword is set
if keyword_set(rep) then var=replicate(var,rep)

; go to the correct location in the file
point_lun,lun,pt

; read the variable
readu,lun,var

; convert the 'char' data type from normal unsigned 'byte' into
; a signed integer
;if strcmp(tp,'char',/fold_case) then begin
;  var = fix(var)
;  idx = where(var gt 127,count)
;  if count ne 0 then var[idx]=var[idx]-256
;endif

end


; function that read the SPE file
;---------------------------------
function read_spe, fname, calib=calib, settings=settings
;+
;  out = read_spe(fname)
;
;  function reads in the SPE-file determined by 'fname' and returns a structure with the spectral data.
;  if the keyword 'calib' is set then the wavelength calibration as defined in '/projects/diagnostics/MAST/MSE/settings'
;  is used rather than the WinSpec wavelength calibration (is more accurate!). In such a case the 'out.wl' field is
;  a [n_pixels x n_channels] array (<=> separate wavelength axis for each channel: taking into account the vertical
;  displacement on the CCD) and the 'out.wlcalib' field will contain the custom calibration coefficient [8 of them and
;  they're not polynomial coefficients!)
;
; :Params:
;    fname   : in, required, type=string
;             filename with if the SPE-file
; :Keywords:
;    calib   : in, type=string
;             perform wavelength calibration as defined in '/projects/diagnostics/MAST/MSE/settings'
;             rather than the WinSpec wavelength calibration (is more accurate!)
;    settings: in, type=string
;             path to the MSE settings directory needed for custom calibration. If not set then following hard coded path is used
;                /projects/diagnostics/MAST/MSE/settings
; :Returns:
;    out     : output structure
;             out.wl   = wavelength vector in nm. [n_pixel x n_channels] array if keyword 'calib' set,
;                        else [n_pixels x 1] array
;             out.t    = time vector (based on readout and exposure time:  t = (1+findgen(n_frames))*(t_read + t_exp)
;             out.ch   = channel vector ( = indgen(n_channels) )
;             out.data = [n_pixels x n_frames x n_channels] array of the calibrated intensity
;             out.wlcalib  = calibration coefficients
;             out.timing.tread  = readout time [s]
;             out.timing.texp   = (extra) exposure time [s]
;             out.timing.ttotal = total exposure time (tread + texp + tclean) [s]
;             out.date = date when data was taken (as 8 digit long integer: yyyymmdd)
;             out.time = local time when data was taken (as 6 digit long integer: hhmmss)
;             out.UTC  = UTC time when data was taken (as 6 digit long integer: hhmmss)
;             out.spect.type    = type of spectrometer (Acton, Spex ,...)
;             out.spect.model   = spectrometer model (e.g. SP275i)
;             out.spect.grating = the type of grating used (number of grooves)
;             out.spect.cwl     = central wavelength the grating is set to [nm]
;             out.spect.slit    = slit width of the entrance slit [um]
;-


  ; open the spe-file for reading
  openr,lun,fname,/get_lun

  ; move through the header and get the relevant info
  ; (the numbers pointing at the relevant location of the data in the file
  ;  can be found in the WinSpec manual, Appendix C. We are using the v2.5-header)

  ; get dimensions of the spectral data
  getd,lun,  42,  'word',n_pixels
  getd,lun,1446,  'long',n_frames
  getd,lun, 656,  'word',n_channels
  ; get the data type of the spectral data
  getd,lun, 108, 'short',datatype
  ; get the exposure and readout times
  getd,lun,  10, 'float',texp
  getd,lun, 672, 'float',tread
  ; get the calibration data
  getd,lun,3101,  'char',polynom_order
  getd,lun,3263,'double',polynom_coeff,rep=6
  if polynom_order lt 6 then wlcalib=polynom_coeff[0:polynom_order] else wlcalib=polynom_coeff
  ; get date/time
  getd,lun,  20,  'char',date,rep=10
  getd,lun, 172,  'char',local,rep=7
  getd,lun, 179,  'char',UTC,rep=7
  ; get spectrometer info
  getd,lun,4043,  'byte',spectype
  getd,lun,4044,  'byte',specmodel
  getd,lun, 650, 'float',grating
  getd,lun,  72, 'float',cwl
  getd,lun, 626, 'float',slit,rep=4


  ; create the wavelength vector
  idx= findgen(n_pixels)+1.
  wl = fltarr(n_pixels)
  for i=0,polynom_order do wl+=polynom_coeff[i] * idx^i

  ; create the timing vector
  tread *= 1e-3          ; convert ms to s 
  ttotal = texp + tread   ; total exposure time
  t  = (1.+findgen(n_frames))*ttotal

  ; create channel vector
  ch = indgen(n_channels)

  ; Convert date and time strings into long integers => easy to sort by date and time
  ; convert date string in 8 digit long integer: yyyymmdd
  day   = 10*(date[0]-48) + (date[1]-48)
  yr    = 1000*(date[5]-48) + 100*(date[6]-48) + 10*(date[7]-48) + (date[8]-48)
  mnstr = strupcase(string(date[2:4]))
  case mnstr of
    'JAN': mn=1
    'FEB': mn=2
    'MAR': mn=3
    'APR': mn=4
    'MAY': mn=5
    'JUN': mn=6
    'JUL': mn=7
    'AUG': mn=8
    'SEP': mn=9
    'OCT': mn=10
    'NOV': mn=11
    'DEC': mn=12
  endcase
  datelng = 10000*long(yr)+100*long(mn)+long(day)
  ; convert times strings in 6 digit long integers: hhmmss
  locallng = 100000*long(local[0]-48) + 10000*long(local[1]-48) $
            + 1000*long(local[2]-48) +   100*long(local[3]-48) $
            +   10*long(local[4]-48) +       long(local[5]-48)
  UTClng   = 100000*long(UTC[0]-48) + 10000*long(UTC[1]-48) $
            + 1000*long(UTC[2]-48) +   100*long(UTC[3]-48) $
            +   10*long(UTC[4]-48) +       long(UTC[5]-48)

  ; Set spectrometer type and model
  if spectype eq 1 then begin
    typestr = 'PI Acton'
    case specmodel of
      1:  modelstr='SP150'
      2:  modelstr='SP275'
      3:  modelstr='SP300'
      4:  modelstr='SP500'
      5:  modelstr='SP500i'
      6:  modelstr='SP750'
      7:  modelstr='SP750i'
      8:  modelstr='AM505'
      9:  modelstr='320PI'
      10: modelstr='InSpectrum'
      11: modelstr='InSight'
    else: modelstr='none'
    endcase
  endif else begin
    typestr = 'SPEX'
    case specmodel of
      1:  modelstr='270M'
      2:  modelstr='ISA TRIAX 180'
      3:  modelstr='ISA TRIAX 190'
      4:  modelstr='ISA TRIAX 320'
      5:  modelstr='232 Retrofit'
      6:  modelstr='ISA TRIAX 550'
      else: modelstr='none'
    endcase
  endelse


  ; read in the actual spectal data:
  case datatype of
  0: data=fltarr(n_pixels,n_channels,n_frames)  ; float datatype
  1: data=lonarr(n_pixels,n_channels,n_frames)  ; long datatype
  2: data=intarr(n_pixels,n_channels,n_frames)  ; integer datatype
  3: data=uintarr(n_pixels,n_channels,n_frames) ; unsigned integer datatype
  endcase
  point_lun,lun,4100  ; got the start of the data
  readu,lun,data      ; and read it in
  dims = size(data,/dim)
  if n_elements(dims) eq 3 then begin
    data = transpose(data,[0,2,1]) ; rearrange the dimensions such that wavelengths are in the first,
                                  ; time is in the second and channels are in the last dimension.
  endif

  ; close the file
  close,lun
  free_lun,lun

  if keyword_set(calib) then begin
    ; wavelength calibration
    ;-----------------------
    ; get settings directory
    if (n_elements(settings) gt 0) && (size(var, /type) eq 7) then begin ; settings keyword exists and is a string
      if file_test(settings[0],/directory) then begin
        settings=settings[0]
      endif else begin
        settings='/projects/diagnostics/MAST/MSE/settings'
      endelse
    endif else begin
      settings='/projects/diagnostics/MAST/MSE/settings'
    endelse
    ; pixel array
    p  = findgen(n_pixels) ; pixels
    p  = -(p-0.5*(n_elements(p)-1))

    ; read in the wavelength calibration file
    wlfile = settings+'/spec_calib/spec_calib.csv'
    print, '* reading wavelength calibration file "'+wlfile+'" ...'
    if ~file_test(wlfile) then begin
      print, '    ERROR: no wavelength calibration file found! No wavelength calibration will be done!'
      goto, NoWL
    endif

    ; wavelength array
    wl = rebin(findgen(n_pixels), n_pixels, n_channels)

    ; get calibration factors
    calib = (read_ascii(wlfile,delim=',')).(0)
    gamma = calib[2, 1]*!dtor
    zeta  = calib[2, 2]*!dtor
    iota  = calib[2, 3]*!dtor
    chi0  = calib[2, 4]*!dtor
    dchi  = calib[2, 5]*!dtor
    w     = calib[2, 6]
    f     = calib[2, 7]
    delta = calib[2, 8]*!dtor
    wlcalib = reform(calib[2,1:8])

    for i=0,n_channels-1 do begin
      chi = chi0 + ch[i]*dchi
      gammas = gamma+zeta
      g      = float(grating)

      ; get theta
      sintheta = - g*cwl*1e-6/2.                             $
                + sin(gammas)                               $
                  * sqrt(2.*(1.+cos(gammas))-(g*cwl*1e-6)^2)$
                  / (2.*(1.+cos(gammas)))
      theta     = asin(sintheta) + iota

      ; get detector part
      tandetect = (p*w*cos(delta)) / (f+p*w*sin(delta))
      detect    = atan(tandetect)

      ; get wl (in [nm])
      wl[*,i] = 1e6*cos(chi)/g * (sin(gamma-(theta+detect)) -  sin(theta))
    endfor
  endif

NoWL:


  ; create data structure
  D = {wl     : wl,                                               $
       t      : t,                                                $
       ch     : ch,                                               $
       data   : data,                                             $
       wlcalib: wlcalib,                                          $
       timing : {tread: tread, texp: texp, ttotal: ttotal},       $
       date   : datelng, time: locallng, UTC: UTClng,             $
       spect  : {type: typestr, model: modelstr, grating: grating,$
                 cwl : cwl,     slit : slit[1]}                   }

  return, D
end
