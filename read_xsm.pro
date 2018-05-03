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

function capof,let,cap
if cap eq 'low' then return,strlowcase(let) else return,strupcase(let)
end

function read_xsm, shot, settings=settings,char=char,cap=cap,dir=dir
if shot le 26251 and shot gt 0 then begin
    default,cap,'high'
    default,char,''
endif else begin
    default,cap,'low'
    default,char,'d'
endelse
;test
;+
;  out = read_xsm(shot)
;
;  reads xsm-data (MSE-spectrometer) and performs wavelength calibration.
;  The returned structure 'out' is (almost) identical to that of the 'read_spe' function.
;  The difference lies in the major radii, fibre bundle no's and NBI-l.o.s. angels (which do
;  not exist in the read_spe return structure), calibration settings (no polynomial), the time (no UTC)
;  and the shot number (which does not exist in the read_spe return structure).
;
; :Params:
;    shot    : in, required, type=long
;             shotnumber for which to read the xsm-data
; :Keywords:
;    settings: in, type=string
;             path to the MSE settings directory. If not set then following hard coded path is used
;                /projects/diagnostics/MAST/MSE/settings
; :Returns:
;    out     : output structure
;             out.wl   = wavelength array [n_pixels x n_channels] in nm (based on wavelength calibration)
;             out.t    = time vector (t[0]=t_trigger+t_exp, t[1]=t[0]+t_read+t_exp, t[2]=t[1]+t_read+t_exp, ...)
;             out.ch   = channel vector
;             out.R    = major radii of the channels (from the MSE patching file)
;             out.md   = fibre bundle numbers of the channels (from the MSE patching file)
;             out.cosb = cosine of the angle between the SS NBI and the line-of-sight.
;                        This is (v.k)/|v| with v the beam velocity vector and k the
;                        emission unit-vector (pointing from the beam towards the collection lens)
;             out.x    = x-coordinates of wheter the lines-of-sight cross the SS NBI
;             out.y    = y-coordinates of wheter the lines-of-sight cross the SS NBI
;             out.xap  = x-coordinate of the aperture of the MSE collection lens
;             out.yap  = y-coordinate of the aperture of the MSE collection lens
;             out.data = [n_pixels x n_frames x n_channels] array of the calibrated intensity
;             out.wlcalib  = coefficients of the calibration [8 x 1 array]
;             out.timing.tread  = readout time [s]
;             out.timing.texp   = (extra) exposure time [s]
;             out.timing.ttotal = total exposure time (tread + texp + tclean) [s]
;             out.date = date when data was taken (as 8 digit long integer: yyyymmdd)
;             out.time = local time when data was taken (as 6 digit long integer: hhmmss)
;             out.shot = shot number
;             out.spect.type    = type of spectrometer (Acton, Spex ,...)
;             out.spect.model   = spectrometer model (e.g. SP275i)
;             out.spect.grating = the type of grating used (number of grooves)
;             out.spect.cwl     = central wavelength the grating is set to [nm]
;             out.spect.slit    = slit width of the entrance slit [um]
;-


  ; initialisation
  ;---------------
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

  ; preparing getdata() for NetCDF => SHOULD BE REMOVED ONCE OFFICIAL IDAM RELEASE CAN HANDLE ATTRIBUTES
  ; NetCDF-4.0 file
  signal_base = 'xsm'
  if shot gt 0 then begin
      ncfile = string(format='("$MAST_DATA/",I03,"/",I0,"/Raw/",A,I06,".nc")',$
                      floor(shot/1000), shot, signal_base, shot)
  endif else begin
      ncfile = string(format='("$MAST_ZSHOT/",A,"z",I05,".nc")',$
                      signal_base,-shot)
  endelse

  if keyword_set(dir) then begin
      ncfile = string(format='(A,"/",A,"z",I05,".nc")',$
                      dir,signal_base,shot)
  endif



  if (file_test(ncfile) eq 0) then begin
    print, 'ERROR: File ' + ncfile + ' does not exist.'
    return, -1
  endif

  ; get the data with getdata()
  source = 'netcdf::'+ncfile


  ; get the raw XSM data
  ;---------------------
  ; shot
  signal  = '/shot'
  data    = getdata(signal,source)
  if data.erc ne 0 then begin
    print, 'ERROR: no data found in file ' + ncfile + '.'
    return, -1
  endif
  shot    = long(data.data[0])
  ; date
  signal  = '/date'
  data    = getdata(signal,source)
  datestr = data.data
  date    = long(strmid(datestr,0,4)+strmid(datestr,5,2)+strmid(datestr,8,2))
  ; time
  signal  = '/time'
  data    = getdata(signal,source)
  timestr = data.data
  time    = long(strmid(timestr,0,2)+strmid(timestr,3,2)+strmid(timestr,6,2))
  ; extra exposure time
  signal  = '/devices/'+char+'2_mse-spec/exposure'
  data    = getdata(signal,source)
  texp    = float(data.data[0])*1e-6
  ; read-out time
  signal  = '/devices/'+char+'2_mse-spec/readout'
  data    = getdata(signal,source)
  tread   = float(data.data[0])*1e-3
  ; total time
  ttotal  = tread+texp
  ; spectrometer type
  signal  = '/devices/'+char+'2_mse-spec/'+capof('S',cap)+'pecType'
  data    = getdata(signal,source)
  type    = string(data.data)
  ; spectrometer model
  signal  = '/devices/'+char+'2_mse-spec/'+capof('S',cap)+'pecModel'
  data    = getdata(signal,source)
  model   = string(data.data)
  ; spectrometer grating
  signal  = '/devices/'+char+'2_mse-spec/'+capof('G',cap)+'ratingDensity'
  data    = getdata(signal,source)
  grating = fix(data.data[0])
  ; CWL [nm]
  signal = '/devices/'+char+'2_mse-spec/wavelength'
  data   = getdata(signal,source)
  cwl    = float(data.data[0])
  ; slit width [micron]
  signal = '/devices/'+char+'2_mse-spec/slit'
  data   = getdata(signal,source)
  slit   = float(data.data[0])

  nch  = 0
  ch   = -1
  data =0.
  for i=0,4 do begin  ; maximum 5 channels for xsm
    signal = string(format='("/",A,"/mse_spec",I1)', signal_base,i)
    tmp    = getdata(signal,source)
    if tmp.erc eq 0 then begin
      if nch eq 0 then begin
        data = tmp.data
        t    = tmp.time
        ch   = [ch,i]
        nch+=1
      endif else begin
        data = [[[data]],[[tmp.data]]]
        ch   = [ch,i]
        nch+=1
      endelse
    endif
  endfor
  if nch eq 0 then begin
    print, 'ERROR: No spectral data found in file ' + ncfile + '.'
    return, -1
  endif
  ch        = ch[1:nch]
  intensity = data

  ; wavelength calibration
  ;-----------------------
  ; pixel array
  p  = findgen(n_elements(data[*,0,0])) ; pixels
  p  = -(p-0.5*(n_elements(p)-1))
  ; wavelength array
  wl = rebin(findgen(n_elements(p)), n_elements(p), nch)

  ; read in the wavelength calibration file
  wlfile = settings+'/spec_calib/spec_calib.csv'
  print, '* reading wavelength calibration file "'+wlfile+'" ...'
  if ~file_test(wlfile) then begin
    print, '    ERROR: no wavelength calibration file found! No wavelength calibration will be done!'
    goto, NoWL
  endif

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

  for i=0,nch-1 do begin
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


NoWL:
  ; Get fibre bundle no's, major radii and cosines with the SS NBI
  ;---------------------------------------------------------------
  md   = fltarr(n_elements(ch))-1.
  R    = fltarr(n_elements(ch))-1.
  cosb = fltarr(n_elements(ch))-1.

  ; read in the patching file
  pfile = settings+'/patching.csv'
  print, '* reading fibre patching from MSE patching file "'+pfile+'" ...'
  if ~file_test(pfile) then begin
    print, '    ERROR: no MSE patching file found! No geometry information will be available!'
    goto, NoR
  endif

  ; reading the patching file:
  ; - patching file structure: 1st line  = header, 2nd to ... lines = data records
  ;                            1st field = shot list, 2nd field = patching beam voltage, 3rd field = radial calibration,
  ;                            4th field = filter measurements, 5th field = spectrometer slit width,
  ;                            6th field = spectrometer central wavelength, 7th field = spectrometer exposure time,
  ;                            8th field = comment field, 9th to ... fields = fibre bundles
  ;                            values in 9th to ... fields = channel/filterscope no. patched to that fibre bundle
  ;                            no value <=> use the value of the previous record
  openr,lun,pfile,/get_lun
  linh = ''
  readf,lun,linh              ; read the header (first line)
  hd   = strsplit(linh,',',/extr,/preserve_null)
  nf   = n_elements(hd)       ; number of fields
  i    = 0
  while not eof(lun) do begin ; loop through the data lines
    ; read the data line
    lin = ''
    readf,lun,lin
    rec = strtrim(strsplit(lin,',',/extr,/preserve_null),2)
    ; initialise the data array on the first line, append the other lines
    if i eq 0 then data = rec else data=[[data],[rec]]
    ; check whether the current record has correct number of fields
    if n_elements(rec) ne nf then begin
      print, "    ERROR: number of fields in the file's header does not match the number of fields in the data record!"
      goto, noR
    endif
    ; if a field of the current record is empty, then copy the data of the same field on the previous record
    ; (ps: NOT for the 'comment' field (j==7)
    for j=0,nf-1 do $
      if (j ne 7) && (strlen(data[j,i]) eq 0) && (i gt 0) then data[j,i] = data[j,i-1]
    ; increment record counter
    i+=1
  endwhile
  nr = i  ; number of records 
  close,lun
  free_lun,lun

  ; compare the desired shot with the shotlist (data[0,*]) to get the correct record
  shl = long(data[0,*])
  idx = where( (shl le shot),count)
  if count eq 0 then begin
    print,shot,format='("    ERROR: shotnumber ",I5," not found in patching file!")
    goto, noR
  endif
  ; the relevant record is the one with the highest shotnumber that is lower or equal to the desired shotnumber.
  idx = idx(n_elements(idx)-1)

  ; the folder with the major radii is given by the 3rd field:
  rposdir = data[2,idx]

  ; spectrometer slit width, wavelength and exposure time
  ;------------------------------------------------------
  slit = uint_extract(data[4,idx],/long)       ; some form of error catching: uint_extract() returns -1 when the string
  if slit ne -1 then slit=float(data[4,idx])$  ; contains no numbers, whereas float() will crash the program
  else slit=0.
  ; spectrometer cwl
  cwl0 = uint_extract(data[5,idx],/long)       ; some form of error catching: uint_extract() returns -1 when the string
  if cwl0 ne -1 then cwl0=float(data[5,idx])$  ; contains no numbers, whereas float() will crash the program
  else cwl0=0.
  ; spectrometer exposure time
  texp0 = uint_extract(data[6,idx],/long)      ; some form of error catching: uint_extract() returns -1 when the string
  if texp0 ne -1 then texp0=float(data[6,idx])$; contains no numbers, whereas float() will crash the program
  else texp0=0.

  if abs(cwl-cwl0) gt 0.01 then begin
    print,format='("    WARNING: the central wavelength that was recorded in the patching file (",F5.1,"nm)")', cwl0
    print,format='("             differs for the actual spectrometer central wavelength (",F5.1,"nm)!")', cwl
    print,format='("             Hence, the patching file log is probably wrong!")'
  endif
  if abs(texp*1e3-texp0) gt 0.01 then begin
    print,format='("    WARNING: the exposure time that was recorded in the patching file (",F5.1,"ms)")', texp0
    print,format='("             differs for the actual spectrometer exposure time (",F5.1,"ms)!")', texp*1e3
    print,format='("             Hence, the patching file log is probably wrong!")'
  endif




  ; all MDs and channels
  mda = fix(strmid(hd[8:nf-1],2))
  cha = data[8:nf-1,idx]          ; includes spectrometer channels (s0 ...), filterscope channels (0 ...) and unpatched (-1)

  ; find the fibre bundles that are connected to the MSE spectrometer
  stag = strmid(cha,0,1) eq 's'
  ns   = fix(total(stag eq 1))
  if ns gt 0 then begin
    chs= fix(strmid(cha[where(stag eq 1)],1))  ; list of patched channels
    mds = mda[where(stag eq 1)]                ; list of fibre bundles connected to the spectrometer
    sidx = sort(chs)
    chs  = chs[sidx]
    mds  = mds[sidx]
    for i=0,nch-1 do begin
      idx=where(ch[i] eq chs, cnt)
      if cnt ne 0 then md[i]=mds[idx[0]]
    endfor
  endif else goto, noR


  ; read in the radial positions file
  ;----------------------------------
  rfile = settings+'/'+rposdir+'/Rpos.csv'
  print, '* reading geometry from MSE geometry file "'+rfile+'" ...'
  if ~file_test(rfile) then begin
    print, '    ERROR: no MSE geometry file found! No geometry information will be available!'
    goto, NoR
  endif
  openr,lun,rfile,/get_lun
  i=0
  ; read (and ignore) the first line (header)
  lin = ''
  readf,lun,lin
  while not eof(lun) do begin ; loop through the data lines
    ; read the data line
    lin = ''
    readf,lun,lin
    rec = strsplit(lin,',',/extr,/preserve_null)
    ; initialise the data array on the first line, append the other lines
    if i eq 0 then data = rec else data=[[data],[rec]]
    i+=1
  endwhile
  nr = i  ; number of records 
  close,lun
  free_lun,lun
  ; convert to floats
  data=float(data)

  ; find fibre bundles of interest
  for i=0,n_elements(md)-1 do begin
    idx = where(data[0,*] eq md[i], cnt)
    if cnt eq 0 then begin
      if md[i] ne -1 then print,format='("    ERROR: fibre bundle MD",I3," was not found in geometry file!")', md[i]
    endif else begin
      R[i]    = data[2,idx]
      cosb[i] = data[3,idx]
    endelse
  endfor


NoR:
  ; get the xy-position of the MSE aperture and the xy-position on the SS NBI
  ;--------------------------------------------------------------------------
  xap = -0.948631
  yap = -2.22780
  x = fltarr(nch)
  y = fltarr(nch)
  xb0 = 0.539
  yb0 =-1.926
  xi  = 85.16*!dtor
  for i=0,nch-1 do begin
    if md[i] ne -1 then begin
      gam = xi + acos(cosb[i])- !pi
      x[i] = (yb0-yap+xap*tan(gam)-xb0*tan(xi))/(tan(gam)-tan(xi))
      y[i] = yap+tan(gam)*(x[i]-xap)
    endif
  endfor

  ; create data structure
  ;----------------------
  out = {wl     : wl,                                         $
         t      : t,                                          $
         ch     : ch,                                         $
         md     : md,                                         $
         R      : R,                                          $
         cosb   : cosb,                                       $
         x      : x,                                          $
         y      : y,                                          $
         xap    : xap,                                        $
         yap    : yap,                                        $
         data   : intensity,                                  $
         wlcalib: wlcalib,                                    $
         timing : {tread: tread, texp: texp, ttotal: ttotal}, $
         date   : date, time: time, shot: shot,               $
         spect  : {type: type, model: model, grating: grating,$
                   cwl : cwl,slit : slit}                     }

  save,out,filename='wavelength_calibration_28101.sav'
  return, out
end
