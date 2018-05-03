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

pro spectrum_compare2, settingfile, inputfiles, labels, prefix, pictdir

; routine overlays the plots of the data in the inputfiles (xdr-format)
; this data can be:
;   - the spatial resolution
;   - the spectra with the different polarised and unpolarised intensities
;   - the polasitation angle
;   - the polarised fraction/intensity
;
; what data is plotted and how it is plotted depends on the settings in the setting file
; (REMARK: some not all settings are ignored: e.g. geometry plot, emission plots)
;
; The 'prefix' is added in front of the title of each plot and (if a Postscript file
; is created) in front of the figures-filename.
;
; Postscript files are saved in the 'pictdir'-directory
;
; v1.0, mdebock 26/07/2007
; v1.1, mdebock 25/07/2007:  * polarised intensity is now: pol. frac * sqrt(total int.)
;                              because that gives a figure-of-merit for the S/N ratio
;                            * some output data is printed to the screen
;
; v1.2, mdebock 31/07/2007:  * Settings for character size, line thickness, ... now depends
;                              on whether we plot to the screen or an EPS-file is written.
;
; v2.0, mdebock 07/08/2008:  * Update to use Stokes vectors as input
;                            * also plots intensity, polarisation angle and polarised fraction profile
;                              for at the CWL (i.e. centre of the sigma or the filter CWL)



;**********************************************************************************
;* READ IN PLOT SETTINGS AND DATA                                                 *
;**********************************************************************************

;----------------------------------------------------------------------------------
; print some info on what's happening to the screen
;----------------------------------------------------------------------------------
; the filenames of the inputfiles without path
inputtxt1 = inputfiles[0]
lastdir   = strpos(inputtxt1, '/', /REVERSE_SEARCH)
inputtxt1 = strmid(inputtxt1,lastdir+1)
inputtxt2 = inputfiles[0]
lastdir   = strpos(inputtxt2, '/', /REVERSE_SEARCH)
inputtxt2 = strmid(inputtxt2,lastdir+1)

print,FORMAT='("* Comparing 2 sets of Stark spectra.")'

;----------------------------------------------------------------------------------
; Read the plot settings
;----------------------------------------------------------------------------------
input = read_setting(settingfile)

;----------------------------------------------------------------------------------
; Read the data
;----------------------------------------------------------------------------------
print,FORMAT='($,"  - reading data from <",A,"> ....")',inputtxt1
restore, inputfiles[0]
chanID1    = chanID
R1         = R
psi1       = psi
R_emis1    = R_emis
psi_emis1  = psi_emis
R_res1     = R_res
psi_res1   = psi_res
lambda1    = lambda
pstokes1   = pstokes
sstokes1   = sstokes
stokes1    = stokes
cwlstokes1 = cwlstokes
print,FORMAT='("... done!")'

; get some inportant data
nchan1   = n_elements(chanID1)
nlambda1 = n_elements(lambda1)
dlambda1 = lambda1[1] - lambda1[0]


print,FORMAT='($,"  - reading data from <",A,"> ....")',inputtxt2
restore, inputfiles[1]
chanID2    = chanID
R2         = R
psi2       = psi
R_emis2    = R_emis
psi_emis2  = psi_emis
R_res2     = R_res
psi_res2   = psi_res
lambda2    = lambda
pstokes2   = pstokes
sstokes2   = sstokes
stokes2    = stokes
cwlstokes2 = cwlstokes
print,FORMAT='("... done!")'


; get some inportant data
nchan2   = n_elements(chanID2[0,*])
nlambda2 = n_elements(lambda1)
dlambda2 = lambda1[1] - lambda1[0]

; the number of channels should be equal
if ~array_equal(chanID1,chanID2) ne 0 then begin
  print, 'ERROR: Both input files need to have the same channels!'
  return
endif
nchan=nchan1

;**********************************************************************************
;* PREPARE FOR PLOTTING OR PRINTING                                               *
;**********************************************************************************

case input.general.plot_print of
  0: return			; no plotting or printing required
  1: begin
       set_plot, 'x'				; plot to the screen
       pseudocol,/on,/quiet			; use pseudocolors
       !p.background=9				; use a white background for all plots
       windowcnt=0
     end
  2: begin
       set_plot, 'x'				; needed to turn on pseudocolors
       pseudocol,/on,/quiet			; use pseudocolors
       set_plot,'PS'				; print to postscript files
       device, scale_factor=1.0, /color,$	; set the postscript device settings
               bits=8, /encapsulated
       ; We don't want white spaces in filenames, so we replace all ' ' by '_' in the
       ; prefix for the ps-filename
       i = 0
       psprefix = prefix
       while 1 do begin
         i = strpos(psprefix,' ',i+1)
         if i eq -1 then break
         strput, psprefix, '_', i
       endwhile
       ; The sizes, line thicknesses, character sizes and thicknesses that produce a nice
       ; plot on the screen are different from that needed to produce a nice eps-file.
       ; These are the conversion parameters: 
       psf  = 0.028			; postscript scale factor
       csf  = 0.85			; character size factor
       thf  = 2.0			; character and line thickness factor
     end
endcase



;**********************************************************************************
;* PLOT THE TOTAL INTENSITY AND THE SPATIAL RESOLUTION                            *
;**********************************************************************************

;----------------------------------------------------------------------------------
; plot the total emission intensity as function of psi and R if requested
;----------------------------------------------------------------------------------
if input.emisres.temisplot  then begin
  ; make the figures title
  ttl = prefix+' - Total emission intensity as function of '+ textoidl('\psi') + ' and R'

  ; set size, character size, character thickness and line thickness
  xsize = 675
  ysize = 875
  csz   = 1.5
  cth   = 1.5
  lth   = 2.0

  ; open a device for plotting the MSE geometry
  if (input.general.plot_print eq 1) then begin
    window, windowcnt, xsize=xsize, ysize=ysize, title=ttl
    windowcnt++
  endif
  if (input.general.plot_print eq 2) then begin	
    device,filename=pictdir+'/'+psprefix+'_total_emission_psi_R.eps',$
           xsize=psf*xsize,ysize=psf*ysize
    csz = csf*csz
    cth = thf*cth
    lth = thf*lth
  endif
  !p.multi = [0,1,2]

  ; calculate the total emission intensity
  dR1       = R1[1]   - R1[0]
  dpsi1     = psi1[1] - psi1[0]
  R_temis1   = total(R_emis1,1)*dR1
  psi_temis1 = total(psi_emis1,1)*dpsi1

  dR2       = R2[1]   - R2[0]
  dpsi2     = psi2[1] - psi2[0]
  R_temis2   = total(R_emis2,1)*dR2
  psi_temis2 = total(psi_emis2,1)*dpsi2

  ; plot the resolution as function of psi
  ttl = prefix+' - Total emission intensity as function of ' + textoidl('\psi')
  xmin = 0
  xmax = 1
  ymin = 0
  ymax = 1.2*max([psi_temis1,psi_temis2])
  plot, [xmin,xmax], [ymin,ymax], /nodata,$
        xs=1, ys=1, color=0,$
        xtitle='normalised flux coordinate '+ textoidl('\psi') +' at Centre of Mass',$
        ytitle='intensity (photons/s)', title=ttl, charsize=csz, charthick=cth
  oplot, psi_res1[2,*], psi_temis1, thick=lth, linestyle=0, psym=-1, symsize=csz, color=3
  oplot, psi_res2[2,*], psi_temis2, thick=lth, linestyle=2, psym=-7, symsize=csz, color=1
  ; indicate what is what
  xyouts, xmax-0.35*(xmax-xmin), ymin+0.9*(ymax-ymin), charsize=csz, charthick=cth, color=3, labels[0]
  xyouts, xmax-0.35*(xmax-xmin), ymin+0.8*(ymax-ymin), charsize=csz, charthick=cth, color=1, labels[1]


  ; plot the resolution as function of R
  ttl = prefix+' - Total emission intensity as function of R'
  xmin = min([R1,R2])
  xmax = max([R1,R2])
  ymin = 0
  ymax = 1.2*max([R_temis1,R_temis2])
  plot, [xmin,xmax], [ymin,ymax], /nodata,$
        xs=1, ys=1, color=0,$
        xtitle='R - Centre of Mass (m)',$
        ytitle='intensity (photons/s)', title=ttl, charsize=csz, charthick=cth
  oplot, R_res1[2,*], R_temis1, thick=lth, linestyle=0, psym=-1, symsize=csz, color=3
  oplot, R_res2[2,*], R_temis2, thick=lth, linestyle=2, psym=-7, symsize=csz, color=1
  ; indicate what is what
  xyouts, xmax-0.35*(xmax-xmin), ymin+0.9*(ymax-ymin), charsize=csz, charthick=cth, color=3, labels[0]
  xyouts, xmax-0.35*(xmax-xmin), ymin+0.8*(ymax-ymin), charsize=csz, charthick=cth, color=1, labels[1]

  ; close the device
  if (input.general.plot_print eq 2) then begin	
    device,/close_file
  endif
  !p.multi = 0

endif
;----------------------------------------------------------------------------------
; plot the spatial resolution as function of psi and R if requested
;----------------------------------------------------------------------------------
if input.emisres.resplot  then begin

  ; make the figures title
  ttl = prefix+' - Spatial resolution as function of '+ textoidl('\psi') + ' and R'

  ; set size, character size, character thickness and line thickness
  xsize = 675
  ysize = 875
  csz   = 1.5
  cth   = 1.5
  lth   = 2.0

  ; open a device for plotting the MSE geometry
  if (input.general.plot_print eq 1) then begin
    window, windowcnt, xsize=xsize, ysize=ysize, title=ttl
    windowcnt++
  endif
  if (input.general.plot_print eq 2) then begin	
    device,filename=pictdir+'/'+psprefix+'_spat_resolution_psi_R.eps',$
           xsize=psf*xsize,ysize=psf*ysize
    csz = csf*csz
    cth = thf*cth
    lth = thf*lth
  endif
  !p.multi = [0,1,2]

  ; plot the resolution as function of psi
  ttl = prefix+' - Spatial resolution as function of ' + textoidl('\psi')
  yttl= 'Spatial resolution ('+textoidl('\psi')+ ') - ' + string(format='(I2)',round(100*psi_res[6,0])) +'% emitted'
  xmin = 0
  xmax = 1
  ymin = 0
  ymax = 1.5*max([psi_res1[5,*],psi_res2[5,*]])
  plot, [xmin,xmax], [ymin,ymax], /nodata,$
        xs=1, ys=1, color=0,$
        xtitle='normalised flux coordinate '+ textoidl('\psi') +' at Centre of Mass',$
        ytitle=yttl, title=ttl, charsize=csz, charthick=cth
  oplot, psi_res1[2,*], psi_res1[5,*], thick=lth, linestyle=0, psym=-1, symsize=csz, color=3
  oplot, psi_res2[2,*], psi_res2[5,*], thick=lth, linestyle=2, psym=-7, symsize=csz, color=1
  ; indicate what is what
  xyouts, xmax-0.35*(xmax-xmin), ymin+0.9*(ymax-ymin), charsize=csz, charthick=cth, color=3, labels[0]
  xyouts, xmax-0.35*(xmax-xmin), ymin+0.8*(ymax-ymin), charsize=csz, charthick=cth, color=1, labels[1]


  ; plot the resolution as function of R
  ttl = prefix+' - Spatial resolution as function of R'
  yttl= 'Spatial resolution (m) - ' + string(format='(I2)',round(100*psi_res[6,0])) +'% emitted'
  xmin = min([R1,R2])
  xmax = max([R1,R2])
  ymin = 0
  ymax = 1.5*max([R_res1[5,*],R_res2[5,*]])
  plot, [xmin,xmax], [ymin,ymax], /nodata,$
        xs=1, ys=1, color=0,$
        xtitle='R - Centre of Mass (m)',$
        ytitle=yttl, title =ttl, charsize=csz, charthick=cth
  oplot, R_res1[2,*], R_res1[5,*], thick=lth, linestyle=0, psym=-1, symsize=csz, color=3
  oplot, R_res2[2,*], R_res2[5,*], thick=lth, linestyle=2, psym=-7, symsize=csz, color=1
  ; indicate what is what
  ; indicate what is what
  xyouts, xmax-0.35*(xmax-xmin), ymin+0.9*(ymax-ymin), charsize=csz, charthick=cth, color=3, labels[0]
  xyouts, xmax-0.35*(xmax-xmin), ymin+0.8*(ymax-ymin), charsize=csz, charthick=cth, color=1, labels[1]

  ; close the device
  if (input.general.plot_print eq 2) then begin	
    device,/close_file
  endif
  !p.multi = 0
endif


;**********************************************************************************
;* PLOT THE SPECTRA                                                               *
;**********************************************************************************

;----------------------------------------------------------------------------------
; Only plot spectral data if requested
;----------------------------------------------------------------------------------
if input.spec.spectrum || input.spec.polangle || input.spec.polfrac then begin

  ; loop through the channels
  for k=0,nchan-1 do begin
    ;----------------------------------------------------------------------------------
    ; If multiplot is selected and some spectral plotting is requested,
    ; then we open a device for plotting here.
    ;----------------------------------------------------------------------------------
    if (input.spec.multiplot) then begin
      nplots=0
      if input.spec.spectrum then nplots++
      if input.spec.polangle then nplots++
      if input.spec.polfrac then nplots++

      ; make the figures title
      ttl = string(format='(" - Channel ",A," at R=",(F4.2),"m / '+textoidl('\psi')+'=",(F4.2))',$
                   chanID1[k],R_res1[2,k],psi_res1[2,k])
      ttl = prefix+ttl

      ; set size, character size, character thickness and line thickness
      xsize = 675.0
      ysize = nplots*285.0
      csz   = 2.5
      csl   = 1.2
      cth   = 1.5
      lth   = 2.0

      ; open a device for plotting the MSE geometry
      if (input.general.plot_print eq 1) then begin
        window, windowcnt, xsize=xsize, ysize=ysize, title=ttl
        windowcnt++
      endif
      if (input.general.plot_print eq 2) then begin
        fname = string(format='("_chan",A,"_R",(F4.2),"_psi",(F4.2),".eps")',$
                      chanID1[k],R_res1[2,k],psi_res1[2,k])
        fname = pictdir+'/'+psprefix+fname
        device,filename=fname, xsize=psf*xsize,ysize=psf*ysize
        csz = csf*csz
        csl = csf*csl
        cth = thf*cth
        lth = thf*lth
      endif
      !p.multi = [0,1,nplots]
    endif

    ;----------------------------------------------------------------------------------
    ; Convert the Stokes vectors in  intensity, polarised fraction,
    ; polarisation angle and norm. S/N
    ;----------------------------------------------------------------------------------
    ; For the first set of data
    nonzero1  = where(stokes1[*,0,k] gt  0.0,count)   ; where the total intensity isn't zero

    tspec1  = fltarr(nlambda1)
    tpspec1 = fltarr(nlambda1)
    ppspec1 = fltarr(nlambda1)
    puspec1 = fltarr(nlambda1)
    spspec1 = fltarr(nlambda1)
    suspec1 = fltarr(nlambda1)
    tpfrac1 = fltarr(nlambda1)+1.0
    SN1     = fltarr(nlambda1)
    pgamma1 = fltarr(nlambda1)
    sgamma1 = fltarr(nlambda1)
    gamma1  = fltarr(nlambda1)

    tspec1[nonzero1]  = stokes1[nonzero1,0,k]
    tpspec1[nonzero1] = sqrt(stokes1[nonzero1,1,k]^2+stokes1[nonzero1,2,k]^2)
    ppspec1[nonzero1] = sqrt(pstokes1[nonzero1,1,k]^2+pstokes1[nonzero1,2,k]^2)
    puspec1[nonzero1] = pstokes1[nonzero1,0,k] - sqrt(pstokes1[nonzero1,1,k]^2+pstokes1[nonzero1,2,k]^2)
    spspec1[nonzero1] = sqrt(sstokes1[nonzero1,1,k]^2+sstokes1[nonzero1,2,k]^2)
    suspec1[nonzero1] = sstokes1[nonzero1,0,k] - sqrt(sstokes1[nonzero1,1,k]^2+sstokes1[nonzero1,2,k]^2)
    tpfrac1[nonzero1] = tpspec1[nonzero1]/tspec1[nonzero1]
    SN1[nonzero1]     = tpspec1[nonzero1]/sqrt(tspec1[nonzero1])

    tpnonzero1 = where(tpspec1 gt 0.0)  ; where the total polarised intensity,
    ppnonzero1 = where(ppspec1 gt 0.0)  ; the polarised pi-intensity
    spnonzero1 = where(spspec1 gt 0.0)  ; and the polarised sigma-intensity isn't zero

    pgamma1[ppnonzero1] = 0.5*atan(pstokes1[ppnonzero1,2,k],pstokes1[ppnonzero1,1,k])*!radeg
    sgamma1[spnonzero1] = 0.5*atan(sstokes1[spnonzero1,2,k],sstokes1[spnonzero1,1,k])*!radeg
    gamma1[tpnonzero1]  = 0.5*atan(stokes1[tpnonzero1,2,k],stokes1[tpnonzero1,1,k])*!radeg

    ; For the second set of data
    nonzero2  = where(stokes2[*,0,k] gt  0.0,count)   ; where the total intensity isn't zero

    tspec2  = fltarr(nlambda2)
    tpspec2 = fltarr(nlambda2)
    ppspec2 = fltarr(nlambda2)
    puspec2 = fltarr(nlambda2)
    spspec2 = fltarr(nlambda2)
    suspec2 = fltarr(nlambda2)
    tpfrac2 = fltarr(nlambda2)+1.0
    SN2     = fltarr(nlambda2)
    pgamma2 = fltarr(nlambda2)
    sgamma2 = fltarr(nlambda2)
    gamma2  = fltarr(nlambda2)

    tspec2[nonzero2]  = stokes2[nonzero2,0,k]
    tpspec2[nonzero2] = sqrt(stokes2[nonzero2,1,k]^2+stokes2[nonzero2,2,k]^2)
    ppspec2[nonzero2] = sqrt(pstokes2[nonzero2,1,k]^2+pstokes2[nonzero2,2,k]^2)
    puspec2[nonzero2] = pstokes2[nonzero2,0,k] - sqrt(pstokes2[nonzero2,1,k]^2+pstokes2[nonzero2,2,k]^2)
    spspec2[nonzero2] = sqrt(sstokes2[nonzero2,1,k]^2+sstokes2[nonzero2,2,k]^2)
    suspec2[nonzero2] = sstokes2[nonzero2,0,k] - sqrt(sstokes2[nonzero2,1,k]^2+sstokes2[nonzero2,2,k]^2)
    tpfrac2[nonzero2] = tpspec2[nonzero2]/tspec2[nonzero2]
    SN2[nonzero2]     = tpspec2[nonzero2]/sqrt(tspec2[nonzero2])

    tpnonzero2 = where(tpspec2 gt 0.0)  ; where the total polarised intensity,
    ppnonzero2 = where(ppspec2 gt 0.0)  ; the polarised pi-intensity
    spnonzero2 = where(spspec2 gt 0.0)  ; and the polarised sigma-intensity isn't zero

    pgamma2[ppnonzero2] = 0.5*atan(pstokes2[ppnonzero2,2,k],pstokes2[ppnonzero2,1,k])*!radeg
    sgamma2[spnonzero2] = 0.5*atan(sstokes2[spnonzero2,2,k],sstokes2[spnonzero2,1,k])*!radeg
    gamma2[tpnonzero2]  = 0.5*atan(stokes2[tpnonzero2,2,k],stokes2[tpnonzero2,1,k])*!radeg

    ;----------------------------------------------------------------------------------
    ; Set the wavelength limits for the plotting
    ;----------------------------------------------------------------------------------
    tmp = max(tspec1,mxidx1)
    if mxidx1 lt 0 then mxidx1=0
    if mxidx1 ge nlambda1 then mxidx1=nlambda1-1
    tmp = max(tspec2,mxidx2)
    if mxidx2 lt 0 then mxidx2=0
    if mxidx2 ge nlambda2 then mxidx2=nlambda2-1

    lambdamin = 0.5*(lambda1[mxidx1]+lambda2[mxidx2])-8.0
    lambdamax = 0.5*(lambda1[mxidx1]+lambda2[mxidx2])+8.0

    ;----------------------------------------------------------------------------------
    ; plot spectrum if requested 
    ;----------------------------------------------------------------------------------
    if input.spec.spectrum then begin

      ;----------------------------------------------------------------------------------
      ; make the figures title
      ;----------------------------------------------------------------------------------
      ttl = string(format='(" - Spectrum ",A," at R=",(F4.2),"m / '+textoidl('\psi')+'=",(F4.2))',$
                  chanID1[k],R_res1[2,k],psi_res1[2,k])
      ttl = prefix+ttl

      ;----------------------------------------------------------------------------------
      ; If multiplot is NOT selected, then we open a device for plotting here.
      ;----------------------------------------------------------------------------------
      if ~(input.spec.multiplot) then begin
        ; set size, character size, character thickness and line thickness
        xsize = 900
        ysize = 600
        csz   = 1.7
        csl   = 1.8
        cth   = 1.5
        lth   = 3.0
        ; open a device for plotting the MSE geometry
        if (input.general.plot_print eq 1) then begin
          window, windowcnt, xsize=xsize, ysize=ysize, title=ttl
          windowcnt++
        endif
        if (input.general.plot_print eq 2) then begin
          fname = string(format='("_spec",A,"_R",(F4.2),"_psi",(F4.2),".eps")',$
                        chanID1[k],R_res1[2,k],psi_res1[2,k])
          fname = pictdir+'/'+psprefix+fname
          device,filename=fname, xsize=psf*xsize,ysize=psf*ysize
          csz = csf*csz
          csl = csf*csl
          cth = thf*cth
          lth = thf*lth
        endif
        !p.multi = 0
      endif


      ;----------------------------------------------------------------------------------
      ; the plotting
      ;----------------------------------------------------------------------------------
      ; set y-limits
      ymin = 0.
      ymax = 1.05 * max([tspec1,tspec2])
      ; axis
      plot, [lambdamin, lambdamax],[ymin,ymax],/nodata,$
            color=0, xs=1, ys=1,$
            xtitle='Wavelength shift (Angstrom)',$
            ytitle='Intensity', $
            title=ttl, charsize=csz, charthick=cth
      ; total spectrum
      oplot,lambda1, tspec1, color=3, thick=lth, linestyle=0
      oplot,lambda2, tspec2, color=1, thick=lth, linestyle=0
      ; pol. pi spectrum
      oplot,lambda1, ppspec1, color=3, thick=lth, linestyle=1
      oplot,lambda2, ppspec2, color=1, thick=lth, linestyle=1
      ; pol. sigma spectrum
      oplot,lambda1, spspec1, color=3, thick=lth, linestyle=2
      oplot,lambda2, spspec2, color=1, thick=lth, linestyle=2
      ; unpol. pi spectrum
      oplot,lambda1, puspec1, color=3, thick=lth, linestyle=3
      oplot,lambda2, puspec2, color=1, thick=lth, linestyle=3
      ; unpol. sigma spectrum
      oplot,lambda1, suspec1, color=3, thick=lth, linestyle=4
      oplot,lambda2, suspec2, color=1, thick=lth, linestyle=4
      ; legend
      xyouts,lambdamin+0.04*(lambdamax-lambdamin),ymin+0.9*(ymax-ymin),$
            charsize=csl, charthick=cth, color=0, '-  total'
      xyouts,lambdamin+0.04*(lambdamax-lambdamin),ymin+0.8*(ymax-ymin),$
            charsize=csl, charthick=cth, color=0, '..   pol. '+textoidl('\pi')
      xyouts,lambdamin+0.04*(lambdamax-lambdamin),ymin+0.7*(ymax-ymin),$
            charsize=csl, charthick=cth, color=0, '-- pol. '+textoidl('\sigma')
      xyouts,lambdamin+0.04*(lambdamax-lambdamin),ymin+0.6*(ymax-ymin),$
            charsize=csl, charthick=cth, color=0, '-.  unpol. '+textoidl('\pi')
      xyouts,lambdamin+0.04*(lambdamax-lambdamin),ymin+0.5*(ymax-ymin),$
            charsize=csl, charthick=cth, color=0, '-..  unpol. '+textoidl('\sigma')
      ; indicate what is what
      xyouts, lambdamax-0.35*(lambdamax-lambdamin),ymin+0.9*(ymax-ymin),$
              charsize=csl, charthick=cth, color=3, labels[0]
      xyouts, lambdamax-0.35*(lambdamax-lambdamin),ymin+0.8*(ymax-ymin),$
              charsize=csl, charthick=cth, color=1, labels[1]

      ;----------------------------------------------------------------------------------
      ; close the device if needed
      ;----------------------------------------------------------------------------------
      if ~(input.spec.multiplot) then begin
        if (input.general.plot_print eq 2) then begin	
          device,/close_file
        endif
        !p.multi = 0
      endif
    endif


    ;----------------------------------------------------------------------------------
    ; plot pol. angle if requested 
    ;----------------------------------------------------------------------------------
    if input.spec.polangle then begin

      ;----------------------------------------------------------------------------------
      ; make the figures title
      ;----------------------------------------------------------------------------------
      ttl = string(format='(" - Polarisation Angle ",A," at R=",(F4.2),"m / '+textoidl('\psi')+'=",(F4.2))',$
                  chanID1[k],R_res1[2,k],psi_res1[2,k])
      ttl = prefix+ttl

      ;----------------------------------------------------------------------------------
      ; If multiplot is NOT selected, then we open a device for plotting here.
      ;----------------------------------------------------------------------------------
      if ~(input.spec.multiplot) then begin
        ; set size, character size, character thickness and line thickness
        xsize = 900
        ysize = 600
        csz   = 1.7
        csl   = 1.8
        cth   = 1.5
        lth   = 3.0
        ; open a device for plotting the MSE geometry
        if (input.general.plot_print eq 1) then begin
          window, windowcnt, xsize=xsize, ysize=ysize, title=ttl
          windowcnt++
        endif
        if (input.general.plot_print eq 2) then begin
          fname = string(format='("_polangle",A,"_R",(F4.2),"_psi",(F4.2),".eps")',$
                        chanID1[k],R_res1[2,k],psi_res1[2,k])
          fname = pictdir+'/'+psprefix+fname
          device,filename=fname, xsize=psf*xsize,ysize=psf*ysize
          csz = csf*csz
          csl = csf*csl
          cth = thf*cth
          lth = thf*lth
        endif
        !p.multi = 0
      endif

      ;----------------------------------------------------------------------------------
      ; the plotting of the pol. angles (only where the polarised intensity is nonzero)
      ;----------------------------------------------------------------------------------
      ; set y-limits
      ymin = -90.0
      ymax = 90.0
      ; axis
      plot, [lambdamin, lambdamax],[ymin,ymax],/nodata,$
            color=0, xs=1, ys=1,$
            xtitle='Wavelength shift (Angstrom)',$
            ytitle='Pol. angle '+textoidl('\gamma')+' (degrees)', $
            title=ttl, charsize=csz, charthick=cth
      ; polarisation angle
      oplot,lambda1[tpnonzero1], gamma1[tpnonzero1], color=3, thick=lth, linestyle=0
      oplot,lambda2[tpnonzero2], gamma2[tpnonzero2], color=1, thick=lth, linestyle=0
      ; indicate what is what
      xyouts, lambdamax-0.35*(lambdamax-lambdamin),ymin+0.9*(ymax-ymin),$
              charsize=csl, charthick=cth, color=3, labels[0]
      xyouts, lambdamax-0.35*(lambdamax-lambdamin),ymin+0.8*(ymax-ymin),$
              charsize=csl, charthick=cth, color=1, labels[1]

      ;----------------------------------------------------------------------------------
      ; close the device if needed
      ;----------------------------------------------------------------------------------
      if ~(input.spec.multiplot) then begin
        if (input.general.plot_print eq 2) then begin	
          device,/close_file
        endif
        !p.multi = 0
      endif
    endif

    ;----------------------------------------------------------------------------------
    ; plot pol. fraction and intensity if requested 
    ;----------------------------------------------------------------------------------
    if input.spec.polfrac then begin

      ;----------------------------------------------------------------------------------
      ; make the figures title
      ;----------------------------------------------------------------------------------
      ttl = string(format='(" - Pol. fraction/SN-ratio ",A," at R=",(F4.2),"m / '+textoidl('\psi')+'=",(F4.2))',$
                  chanID1[k],R_res1[2,k],psi_res1[2,k])
      ttl = prefix+ttl

      ;----------------------------------------------------------------------------------
      ; If multiplot is NOT selected, then we open a device for plotting here.
      ;----------------------------------------------------------------------------------
      if ~(input.spec.multiplot) then begin
        ; set size, character size, character thickness and line thickness
        xsize = 900
        ysize = 600
        csz   = 1.7
        csl   = 1.8
        cth   = 1.5
        lth   = 3.0
        ; open a device for plotting the MSE geometry
        if (input.general.plot_print eq 1) then begin
          window, windowcnt, xsize=xsize, ysize=ysize, title=ttl
          windowcnt++
        endif
        if (input.general.plot_print eq 2) then begin
          fname = string(format='("_polfrac",A,"_R",(F4.2),"_psi",(F4.2),".eps")',$
                        chanID1[k],R_res1[2,k],psi_res1[2,k])
          fname = pictdir+'/'+psprefix+fname
          device,filename=fname, xsize=psf*xsize,ysize=psf*ysize
          csz = csf*csz
          csl = csf*csl
          cth = thf*cth
          lth = thf*lth
        endif
        !p.multi = 0
      endif


      ;----------------------------------------------------------------------------------
      ; the plotting of the pol. fraction (only where the polarised intensity is nonzero)
      ;----------------------------------------------------------------------------------
      ; set y-limits
      ymin = 0.0
      ymax = 1.0
      ; axis
      plot, [lambdamin, lambdamax],[ymin,ymax],/nodata,$
            color=0, xs=1, ys=1,$
            xtitle='Wavelength shift (Angstrom)',$
            ytitle='Polarised fraction/norm. S/N-ratio', $
            title=ttl, charsize=csz, charthick=cth
      ; polarised fraction
      oplot, lambda1[tpnonzero1], tpfrac1[tpnonzero1], color=3, linestyle=0, thick=lth
      oplot, lambda2[tpnonzero2], tpfrac2[tpnonzero2], color=1, linestyle=0, thick=lth
      ; normalised polarised intensity
      maxSN = max([max(SN1),max(SN2)])
      oplot, lambda1[tpnonzero1], SN1[tpnonzero1]/maxSN, color=3, linestyle=2, thick=lth
      oplot, lambda2[tpnonzero2], SN2[tpnonzero2]/maxSN, color=1, linestyle=2, thick=lth
      ; legend
      xyouts,lambdamin+0.04*(lambdamax-lambdamin),ymin+0.9*(ymax-ymin),$
            charsize=csl, charthick=cth, color=0, '-   pol. fraction'
      xyouts,lambdamin+0.04*(lambdamax-lambdamin),ymin+0.8*(ymax-ymin),$
            charsize=csl, charthick=cth, color=0, '-- norm. S/N-ratio'
      ; indicate what is what
      xyouts, lambdamax-0.35*(lambdamax-lambdamin),ymin+0.9*(ymax-ymin),$
              charsize=csl, charthick=cth, color=3, labels[0]
      xyouts, lambdamax-0.35*(lambdamax-lambdamin),ymin+0.8*(ymax-ymin),$
              charsize=csl, charthick=cth, color=1, labels[1]


      ;----------------------------------------------------------------------------------
      ; close the device if needed
      ;----------------------------------------------------------------------------------
      if ~(input.spec.multiplot) then begin
        if (input.general.plot_print eq 2) then begin	
          device,/close_file
        endif
        !p.multi = 0
      endif
    endif

    ;----------------------------------------------------------------------------------
    ; close the device if needed
    ;----------------------------------------------------------------------------------
    if (input.spec.multiplot) then begin
      if (input.general.plot_print eq 2) then begin	
        device,/close_file
      endif
      !p.multi = 0
    endif

  endfor
endif


;**********************************************************************************
;* PLOT THE PROFILES                                                              *
;**********************************************************************************


;----------------------------------------------------------------------------------
; Only plot spectral data if requested
;----------------------------------------------------------------------------------
if input.profile.intensity || input.profile.polangle || input.profile.polfrac then begin

  ;----------------------------------------------------------------------------------
  ; If multiplot is selected and some spectral plotting is requested,
  ; then we open a device for plotting here.
  ;----------------------------------------------------------------------------------
  if (input.spec.multiplot) then begin
    nplots=0
    if input.spec.spectrum then nplots++
    if input.spec.polangle then nplots++
    if input.spec.polfrac then nplots++

    ; make the figures title
    ttl = ' - Profiles at CWL'
    ttl = prefix+ttl

    ; set size, character size, character thickness and line thickness
    xsize = 675.0
    ysize = nplots*285.0
    csz   = 2.5
    csl   = 1.2
    cth   = 1.5
    lth   = 2.0

    ; open a device for plotting the spectra
    if (input.general.plot_print eq 1) then begin
      window, windowcnt, xsize=xsize, ysize=ysize, title=ttl
      windowcnt++
    endif
    if (input.general.plot_print eq 2) then begin
      fname = '_profilesCWL.eps'
      fname = pictdir+'/'+psprefix+fname
      device,filename=fname, xsize=psf*xsize,ysize=psf*ysize
      csz = csf*csz
      csl = csf*csl
      cth = thf*cth
      lth = thf*lth
    endif
    !p.multi = [0,1,nplots]
  endif

  ;----------------------------------------------------------------------------------
  ; Convert the Stokes vectors in intensity, polarised fraction,
  ; polarisation angle and norm. S/N 
  ;----------------------------------------------------------------------------------
  ; For the first set of data
  intcwl1   = reform(cwlstokes1[0,*])
  pfraccwl1 = reform(sqrt(cwlstokes1[1,*]^2+cwlstokes1[2,*]^2)/cwlstokes1[0,*])
  SNcwl1    = reform(sqrt(cwlstokes1[1,*]^2+cwlstokes1[2,*]^2/cwlstokes1[0,*]))
  gammacwl1 = 0.5*reform(atan(cwlstokes1[2,*],cwlstokes1[1,*]))*!radeg
  CWL1      = cwlstokes1[3,*]
  ; For the second set of data
  intcwl2   = reform(cwlstokes2[0,*])
  pfraccwl2 = reform(sqrt(cwlstokes2[1,*]^2+cwlstokes2[2,*]^2)/cwlstokes2[0,*])
  SNcwl2    = reform(sqrt(cwlstokes2[1,*]^2+cwlstokes2[2,*]^2/cwlstokes2[0,*]))
  gammacwl2 = 0.5*reform(atan(cwlstokes2[2,*],cwlstokes2[1,*]))*!radeg
  CWL2      = cwlstokes2[3,*]

  ;----------------------------------------------------------------------------------
  ; plot intensity profile if requested 
  ;----------------------------------------------------------------------------------
  if input.profile.intensity then begin

    ;----------------------------------------------------------------------------------
    ; make the figures title
    ;----------------------------------------------------------------------------------
    ttl = ' - Intensity at CWL'
    ttl = prefix+ttl

    ;----------------------------------------------------------------------------------
    ; If multiplot is NOT selected, then we open a device for plotting here.
    ;----------------------------------------------------------------------------------
    if ~(input.profile.multiplot) then begin
      ; set size, character size, character thickness and line thickness
      xsize = 900
      ysize = 600
      csz   = 1.7
      csl   = 1.8
      cth   = 1.5
      lth   = 3.0
      ; open a device for plotting the MSE geometry
      if (input.general.plot_print eq 1) then begin
        window, windowcnt, xsize=xsize, ysize=ysize, title=ttl
        windowcnt++
      endif
      if (input.general.plot_print eq 2) then begin
        fname = '_intensityCWL.eps")'
        fname = pictdir+'/'+psprefix+fname
        device,filename=fname, xsize=psf*xsize,ysize=psf*ysize
        csz = csf*csz
        csl = csf*csl
        cth = thf*cth
        lth = thf*lth
      endif
      !p.multi = 0
    endif

    ;----------------------------------------------------------------------------------
    ; the plotting
    ;----------------------------------------------------------------------------------
    ; set limits
    xmin = min(R)
    xmax = max(R)
    ymin = 0.
    ymax = 1.05*max([max(intcwl1),max(intcwl2)])
    ; axis
    plot, [xmin, xmax],[ymin,ymax],/nodata,$
          color=0, xs=1, ys=1,$
          xtitle='R - Centre of Mass (m)',$
          ytitle='Intensity', $
          title=ttl, charsize=csz, charthick=cth
    ; intensity
    oplot, R_res1[2,*], intcwl1, color=3, linestyle=0, thick=lth, psym=-1
    oplot, R_res2[2,*], intcwl2, color=1, linestyle=0, thick=lth, psym=-7
    ; indicate what is what
    xyouts, xmax-0.35*(xmax-xmin),ymin+0.9*(ymax-ymin),$
            charsize=csl, charthick=cth, color=3, labels[0]
    xyouts, xmax-0.35*(xmax-xmin),ymin+0.8*(ymax-ymin),$
            charsize=csl, charthick=cth, color=1, labels[1]

    ;----------------------------------------------------------------------------------
    ; close the device if needed
    ;----------------------------------------------------------------------------------
    if ~(input.profile.multiplot) then begin
      if (input.general.plot_print eq 2) then begin	
        device,/close_file
      endif
      !p.multi = 0
    endif

  endif


  ;----------------------------------------------------------------------------------
  ; plot pol. angle if requested 
  ;----------------------------------------------------------------------------------
  if input.profile.polangle then begin

    ;----------------------------------------------------------------------------------
    ; make the figures title
    ;----------------------------------------------------------------------------------
    ttl = ' - Polarisation Angle at CWL'
    ttl = prefix+ttl

    ;----------------------------------------------------------------------------------
    ; If multiplot is NOT selected, then we open a device for plotting here.
    ;----------------------------------------------------------------------------------
    if ~(input.profile.multiplot) then begin
      ; set size, character size, character thickness and line thickness
      xsize = 900
      ysize = 600
      csz   = 1.7
      csl   = 1.8
      cth   = 1.5
      lth   = 3.0
      ; open a device for plotting the MSE geometry
      if (input.general.plot_print eq 1) then begin
        window, windowcnt, xsize=xsize, ysize=ysize, title=ttl
        windowcnt++
      endif
      if (input.general.plot_print eq 2) then begin
        fname = '_polangleCWL.eps'
        fname = pictdir+'/'+psprefix+fname
        device,filename=fname, xsize=psf*xsize,ysize=psf*ysize
        csz = csf*csz
        csl = csf*csl
        cth = thf*cth
        lth = thf*lth
      endif
      !p.multi = 0
    endif

    ;----------------------------------------------------------------------------------
    ; the plotting of the pol. angle profile
    ;----------------------------------------------------------------------------------
    ; set limits
    xmin = min(R)
    xmax = max(R)
    ymin = -90.0
    ymax = 90.0
    ; axis
    plot, [xmin, xmax],[ymin,ymax],/nodata,$
          color=0, xs=1, ys=1,$
          xtitle='R - Centre of Mass (m)',$
          ytitle='Pol. angle '+textoidl('\gamma')+' (degrees)', $
          title=ttl, charsize=csz, charthick=cth
    ; pol. angle
    oplot, R_res1[2,*], gammacwl1, color=3, linestyle=0, thick=lth, psym=-1
    oplot, R_res2[2,*], gammacwl2, color=1, linestyle=0, thick=lth, psym=-7
    ; indicate what is what
    xyouts, xmax-0.35*(xmax-xmin),ymin+0.9*(ymax-ymin),$
            charsize=csl, charthick=cth, color=3, labels[0]
    xyouts, xmax-0.35*(xmax-xmin),ymin+0.8*(ymax-ymin),$
            charsize=csl, charthick=cth, color=1, labels[1]

    ;----------------------------------------------------------------------------------
    ; close the device if needed
    ;----------------------------------------------------------------------------------
    if ~(input.profile.multiplot) then begin
      if (input.general.plot_print eq 2) then begin	
        device,/close_file
      endif
      !p.multi = 0
    endif

  endif

  ;----------------------------------------------------------------------------------
  ; plot pol. fraction and norm. S/N ratio if requested 
  ;----------------------------------------------------------------------------------
  if input.profile.polfrac then begin

    ;----------------------------------------------------------------------------------
    ; make the figures title
    ;----------------------------------------------------------------------------------
    ttl = ' - Pol. fraction/SN-ratio at CWL'
    ttl = prefix+ttl

    ;----------------------------------------------------------------------------------
    ; If multiplot is NOT selected, then we open a device for plotting here.
    ;----------------------------------------------------------------------------------
    if ~(input.profile.multiplot) then begin
      ; set size, character size, character thickness and line thickness
      xsize = 900
      ysize = 600
      csz   = 1.7
      csl   = 1.8
      cth   = 1.5
      lth   = 3.0
      ; open a device for plotting the MSE geometry
      if (input.general.plot_print eq 1) then begin
        window, windowcnt, xsize=xsize, ysize=ysize, title=ttl
        windowcnt++
      endif
      if (input.general.plot_print eq 2) then begin
        fname = '_polfracCWL.eps'
        fname = pictdir+'/'+psprefix+fname
        device,filename=fname, xsize=psf*xsize,ysize=psf*ysize
        csz = csf*csz
        csl = csf*csl
        cth = thf*cth
        lth = thf*lth
      endif
      !p.multi = 0
    endif

    ;----------------------------------------------------------------------------------
    ; the plotting of the pol. fraction (only where the polarised intensity is nonzero)
    ;----------------------------------------------------------------------------------
    ; set limits
    xmin = min(R)
    xmax = max(R)
    ymin = 0.0
    ymax = 1.0
    ; axis
    plot, [xmin, xmax],[ymin,ymax],/nodata,$
          color=0, xs=1, ys=1,$
          xtitle='R - Centre of Mass (m)',$
          ytitle='Polarised fraction/ norm. S/N-ratio', $
          title=ttl, charsize=csz, charthick=cth
    ; pol. fraction
    oplot, R_res1[2,*], pfraccwl1, color=3, linestyle=0, thick=lth, psym=-1
    oplot, R_res2[2,*], pfraccwl2, color=1, linestyle=0, thick=lth, psym=-7
    ; normalised S/N ratio
    oplot, R_res1[2,*], SNcwl1/max(SNcwl1), color=3, linestyle=2, thick=lth, psym=-1
    oplot, R_res2[2,*], SNcwl2/max(SNcwl2), color=1, linestyle=2, thick=lth, psym=-7
    ; legend
    xyouts,xmin+0.04*(xmax-xmin),ymin+0.9*(ymax-ymin),$
           charsize=csl, charthick=cth, color=0, '-   pol. fraction'
    xyouts,xmin+0.04*(xmax-xmin),ymin+0.8*(ymax-ymin),$
           charsize=csl, charthick=cth, color=0, '-- norm. S/N ratio'
    ; indicate what is what
    xyouts, xmax-0.35*(xmax-xmin),ymin+0.9*(ymax-ymin),$
            charsize=csl, charthick=cth, color=3, labels[0]
    xyouts, xmax-0.35*(xmax-xmin),ymin+0.8*(ymax-ymin),$
            charsize=csl, charthick=cth, color=1, labels[1]

    ;----------------------------------------------------------------------------------
    ; close the device if needed
    ;----------------------------------------------------------------------------------
    if ~(input.profile.multiplot) then begin
      if (input.general.plot_print eq 2) then begin	
        device,/close_file
      endif
      !p.multi = 0
    endif

  endif

  ;----------------------------------------------------------------------------------
  ; close the device if needed
  ;----------------------------------------------------------------------------------
  if (input.profile.multiplot) then begin
    if (input.general.plot_print eq 2) then begin	
      device,/close_file
    endif
    !p.multi = 0
  endif

endif



;----------------------------------------------------------------------------------
; Reset to default plotting
;----------------------------------------------------------------------------------
if (input.general.plot_print eq 2) then begin	
  device, scale_factor=1.0, color=0,$	; resetting the default postscript device settings
          encapsulated=0
endif
print, ''
set_plot, 'x'	; plot to the screen
!p.multi = 0



end
