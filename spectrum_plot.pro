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

pro spectrum_plot, settingfile, inputfile, prefix, pictdir

; routine plots the data form the given inputfile (xdr-format)
; this data can be:
;   - the 3D geometry of the MSE system
;   - the emission in the (R,Z)-plane, as function of R andof psi
;   - the spatial resolution
;   - the spectra with the different polarised and unpolarised intensities
;   - the polasitation angle
;   - the polarised fraction and intensity
;
; what data is plotted and how it is plotted depends on the settings in the setting file
;
; The 'prefix' is added in front of the title of each plot and (if a Postscript file
; is created) in front of the figures-filename.
;
; Postscript files are saved in the 'pictdir'-directory
;
; v1.0, mdebock 18/07/2007
;
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
; v2.1, mdebock 18/07/2011:  * Better determination of the wavelength range
;                            * Better handling of colours and EPS/X-window plotting with the init_graphic and truecolor functions
; v2.2, daussems 29/09/2011: * Ideal 1D MSE angle now indicated + deviation from full modelling
;       mdebock              * Position optimal pi-red, sigma and pi-blue now determined outside spectrum plot
;                            * Extra profile plots at optimal pi-red, sigma and pi-blue (CWL is of course still there)
;                            * Windows/Unix compatibility ensured
; v2.3, pgeelen 09/02/2012:  * Update to use 4D Stokes vectors as input 
;       mdebock       
;

;**********************************************************************************
;* READ IN PLOT SETTINGS AND DATA                                                 *
;**********************************************************************************
;----------------------------------------------------------------------------------
; Distinguish between Windows and UNIX directory separators
;----------------------------------------------------------------------------------
if strcmp( !version.os,'Win32',/fold_case) then sep='\' else sep='/'


;----------------------------------------------------------------------------------
; we will put the info that is printed to the screen also in a html file 
; (is easy to write and can easily be converted in postscript or pdf)
;----------------------------------------------------------------------------------
; the filename of the inputfile without path or extension
inputtxt = inputfile
print,'input text', inputfile
lastdir  = strpos(inputtxt, sep, /REVERSE_SEARCH)
inputtxt = strmid(inputtxt,lastdir+1)
ext      = strpos(inputtxt, '.', /REVERSE_SEARCH)
inputtxt = strmid(inputtxt,0,ext)
print,'inputtxt',inputtxt

; name of the textfile 
txtfile = pictdir+sep+inputtxt+'.html'
; open the file
get_lun, funit
openw, funit, txtfile
printf, funit, '<html>'
printf, funit, '<head>'
printf, funit, '<title>Output data of '+inputtxt+'.xdr</title>'
printf, funit, '<style type="text/css">'
printf, funit, '<!--'
printf, funit, ' h1 {font-size: 150%; }'
printf, funit, ' h2 {font-size: 125%; }'
printf, funit, '-->'
printf, funit, '</style>'
printf, funit, '</head>'
printf, funit, '<body>'
printf, funit, '<h1>Output data "'+inputtxt+'.xdr"</h1>

;----------------------------------------------------------------------------------
; print some info on what's happening to the screen
;----------------------------------------------------------------------------------
print,FORMAT=  '("* Plotting ",(A-25))', inputtxt+'.xdr'
print,FORMAT='($,"  - reading data ....")'

;----------------------------------------------------------------------------------
; Read the plot settings
;----------------------------------------------------------------------------------
input = read_setting(settingfile)

;----------------------------------------------------------------------------------
; Read the data
;----------------------------------------------------------------------------------
restore, inputfile
print,FORMAT='("... done!")'

; get some inportant data
nchan   = n_elements(gp_xyz[0,0,*])
gp_n    = n_elements(gp_xyz[0,*,0])
nR      = n_elements(R)
nZ      = n_elements(Z)
npsi    = n_elements(psi)
nlambda = n_elements(lambda)
dlambda = lambda[1] - lambda[0]


;**********************************************************************************
;* PREPARE FOR PLOTTING OR PRINTING                                               *
;**********************************************************************************
case input.general.plot_print of
  0: return			; no plotting or printing required
  1: begin
       epsfig = 0
       ; font 'complex roman' the X plots
       fontname = 'complex roman'
     end
  2: begin
       epsfig=1
       ; We don't want white spaces in filenames, so we replace all ' ' by '_' in the
       ; prefix for the ps-filename
       i = 0
       psprefix = prefix
       while 1 do begin
         i = strpos(psprefix,' ',i+1)
         if i eq -1 then break
         strput, psprefix, '_', i
       endwhile
       ; font 'times' the EPS plots
       fontname = 'times'
     end
endcase

;**********************************************************************************
;* PLOT MSE SYSTEM GEOMETRY IN 3D                                                 *
;**********************************************************************************
if input.geom.geomplot then begin

  ;----------------------------------------------------------------------------------
  ; open a device and plot the 3D axis
  ;----------------------------------------------------------------------------------
  ; get centre, zoom and rotation for the 3D axis
  centre   = input.geom.centre
  zoom     = input.geom.zoom
  rotation = input.geom.rotation


  ; set aspect ratio, figure zoom, fontsize and line thickness
  aspect   = 1
  figzoom  = 1.2
  fontsize = 16
  lth      = 2.0
 
  ; make the figures title
  ttl = prefix+' - MSE system geometry in 3D'
  ; open a device for plotting the MSE geometry
  if (input.general.plot_print eq 1) then begin	
    figID = init_graphic(aspect=aspect, fontsize=fontsize, fontname=fontname, zoom=figzoom, title=ttl)
  endif
  if (input.general.plot_print eq 2) then begin	
    fname = pictdir+'/'+psprefix+'_MSE_geometry.eps'
    figID = init_graphic(aspect=aspect, fontsize=fontsize, fontname=fontname, zoom=figzoom, eps=fname, bits_per_pixel=8)
  endif
  !p.multi=0
  !p.background=truecolor('white')	; use a white background
  !p.color     =truecolor('black')	; use a black foreground
  
  ; plot a 3D-axis
  plot_axis, centre, zoom, rotation, ttl, lth

  ;----------------------------------------------------------------------------------
  ; plot a beam wireframe is requested
  ;----------------------------------------------------------------------------------
  if input.geom.beamplot then begin
    ; the beam range is chosen to be 20% more than the distance 
    ; between the first maximum and minimum radius
    dRmax = quadroots([ B_xyz[0]^2 + B_xyz[1]^2 - max(R)^2 ,$
                        2*(B_xyz[0]*B_vec[0] + B_xyz[1]*B_vec[1]),$
                        B_vec[0]^2 + B_vec[1]^2])
    dRmax = min(dRmax)
    dRmin = quadroots([ B_xyz[0]^2 + B_xyz[1]^2 - min(R)^2 ,$
                        2*(B_xyz[0]*B_vec[0] + B_xyz[1]*B_vec[1]),$
                        B_vec[0]^2 + B_vec[1]^2])
    dRmin = min(dRmin)

    B_range = [dRmax - 0.1*(dRmin-dRmax),$
               dRmin + 0.1*(dRmin-dRmax) ]
    plot_beam, B_xyz, B_vec, B_range, B_w
  endif

  ;----------------------------------------------------------------------------------
  ; plot the grid points if requested
  ;----------------------------------------------------------------------------------
  colors = truecolor(/help)
  if input.geom.gpplot then begin
    for k=1,nchan do begin
      plots, gp_xyz[0,*,nchan-k], gp_xyz[1,*,nchan-k], gp_xyz[2,*,nchan-k],$
             color=truecolor(colors[k mod n_elements(colors)]), psym=1, thick=0.5*lth, symsize=csz, /T3D, noclip=0
    endfor
  endif

  ;----------------------------------------------------------------------------------
  ; plot the B-field (T) if requested
  ;----------------------------------------------------------------------------------
  if input.geom.Bfldplot then begin
    for k=1,nchan do begin
      for gp_c=0,gp_n-1 do begin
        plots, [gp_xyz[0,gp_c,nchan-k], gp_xyz[0,gp_c,nchan-k]+gp_Bfld[0,gp_c,nchan-k]],$
               [gp_xyz[1,gp_c,nchan-k], gp_xyz[1,gp_c,nchan-k]+gp_Bfld[1,gp_c,nchan-k]],$
               [gp_xyz[2,gp_c,nchan-k], gp_xyz[2,gp_c,nchan-k]+gp_Bfld[2,gp_c,nchan-k]],$
               color=truecolor('blue'), thick=lth, /T3D, noclip=0
      endfor
    endfor
  endif

  ;----------------------------------------------------------------------------------
  ; plot beam particle velocity (1e7 m/s) if requested
  ;----------------------------------------------------------------------------------
  if input.geom.velplot then begin
    for k=1,nchan do begin
      for gp_c=0,gp_n-1 do begin
        plots, [gp_xyz[0,gp_c,nchan-k], gp_xyz[0,gp_c,nchan-k]+B_v0*gp_vel[0,gp_c,nchan-k]*1e-7],$
               [gp_xyz[1,gp_c,nchan-k], gp_xyz[1,gp_c,nchan-k]+B_v0*gp_vel[1,gp_c,nchan-k]*1e-7],$
               [gp_xyz[2,gp_c,nchan-k], gp_xyz[2,gp_c,nchan-k]+B_v0*gp_vel[2,gp_c,nchan-k]*1e-7],$
               color=truecolor('red'), thick=lth, /T3D, noclip=0
      endfor
    endfor
  endif

  ;----------------------------------------------------------------------------------
  ; plot E-field (1e7 V/m) if requested
  ;----------------------------------------------------------------------------------
  if input.geom.Efldplot then begin
    for k=1,nchan do begin
      for gp_c=0,gp_n-1 do begin
        plots, [gp_xyz[0,gp_c,nchan-k], gp_xyz[0,gp_c,nchan-k]+gp_Efld[0,gp_c,nchan-k]*1e-7],$
               [gp_xyz[1,gp_c,nchan-k], gp_xyz[1,gp_c,nchan-k]+gp_Efld[1,gp_c,nchan-k]*1e-7],$
               [gp_xyz[2,gp_c,nchan-k], gp_xyz[2,gp_c,nchan-k]+gp_Efld[2,gp_c,nchan-k]*1e-7],$
               color=truecolor('green'), thick=lth, /T3D, noclip=0
      endfor
    endfor
  endif


  ;----------------------------------------------------------------------------------
  ; plot beam emission if requested
  ;----------------------------------------------------------------------------------
  if input.geom.emisplot then begin
    ; load a colour tabel
    if !d.name eq 'X' then begin
      device, get_decomp=decomp
      device, decomp=0
    endif
    loadct, 5, /silent
    n_colors = !d.table_size
    ; get the colour scaling
    emismax = max(gp_emis)
    emismin = min(gp_emis)
    if abs(emismax-emismin) gt 1e-6 then cfactor  = (n_colors-1.)/(emismax-emismin)$
    else cfactor = (n_colors-1.)/(2.*emismax)
    for k=1,nchan do begin
      for gp_c=0l,gp_n-1 do begin
        coloridx = round(cfactor*(gp_emis[gp_c,nchan-k]-emismin))
        plots, gp_xyz[0,gp_c,nchan-k],gp_xyz[1,gp_c,nchan-k],gp_xyz[2,gp_c,nchan-k],$
               color=coloridx, psym=1, thick=0.5*lth, symsize=csz, /T3D
      endfor
    endfor
    if !d.name eq 'X' then device, decomp=decomp
  endif

  ;----------------------------------------------------------------------------------
  ; close the figure
  ;----------------------------------------------------------------------------------
  tmp = reset_graphic(figID, eps=epsfig)
  !p.multi = 0

endif

;**********************************************************************************
;* PLOT THE EMISSION INTENSITY AND SPATIAL RESOLUTION                             *
;**********************************************************************************

;----------------------------------------------------------------------------------
; plot the beam emission in a RZ-plane if requested
;----------------------------------------------------------------------------------
if input.emisres.emisRZplot then begin

  ; make the figures title
  ttl = prefix+' - Beam emission in the RZ-plane'

  ; set aspect ratio, figure zoom, fontsize and line thickness
  if nR ge nZ then aspect   = 1.6 else aspect = 0.7
  figzoom  = 2.0
  fontsize = 14
  lth      = 2.0

  ; open a device for plotting the RZ emission
  if (input.general.plot_print eq 1) then begin	
    figID = init_graphic(aspect=aspect, fontsize=fontsize, fontname=fontname, zoom=figzoom, title=ttl)
    !p.background=truecolor('white')	; use a white background
    !p.color     =truecolor('black')	; use a black foreground
  endif
  if (input.general.plot_print eq 2) then begin	
    fname = pictdir+'/'+psprefix+'_RZ_emission.eps'
    figID = init_graphic(aspect=aspect, fontsize=fontsize, fontname=fontname, /PSfont, zoom=figzoom, eps=fname, bits_per_pixel=8)
  endif
  !p.multi=0
  !p.background=truecolor('white')	; use a white background
  !p.color     =truecolor('black')	; use a black foreground
 
  ; add up the RZ-emission data for all each channels
  ; the colorscale
  RZ_total =fltarr(nR,nZ)
  for k=0,nchan-1 do begin
    if max(RZ_emis[*,*,k]) gt 0.0 then begin
      RZ_total += RZ_emis[*,*,k]
    endif
  endfor

  ; plot the shaded contour plot
  cp_shaded, RZ_total, R, Z, /isotropic, xtitle='R [m]', ytitle='Z [m]', ztitle='Intensity [photons/s/m^2]',$
            title=ttl, zr=[0,max(RZ_total)],ctable=3, /invert, /showscale ; pixels=4096

  ; overlay the fluxsurfaces and plot the magnetic axis
  contour,RZpsi,R,Z, levels=0.1*(1+findgen(10)),c_labels=1+intarr(10),$
          /overplot,color=0,thick=lth,xs=1,ys=1
  plots,Rm,0.0, psym=1, color=0,thick=lth

  ; close the figure
  tmp = reset_graphic(figID, eps=epsfig)
  !p.multi = 0
endif 

;----------------------------------------------------------------------------------
; plot the beam emission as function of psi and R if requested
;----------------------------------------------------------------------------------
if input.emisres.emisPsiRplot  then begin
  ; make the figures title
  ttl = prefix+' - Beam emission intensity as function of ' + textoidl('\psi') + ' and R'

  ; set aspect ratio, figure zoom, fontsize and line thickness
  aspect   = 0.8
  figzoom  = 1.2
  fontsize = 12
  lth      = 2.0
 
  ; open a device for plotting the emission as function of R and psi
  if (input.general.plot_print eq 1) then begin	
    figID = init_graphic(aspect=aspect, fontsize=fontsize, fontname=fontname, zoom=figzoom, title=ttl)
  endif
  if (input.general.plot_print eq 2) then begin	
    fname = pictdir+'/'+psprefix+'_emission_psi_R.eps'
    figID = init_graphic(aspect=aspect, fontsize=fontsize, fontname=fontname, /PSfont, zoom=figzoom, eps=fname, bits_per_pixel=8)
  endif
  !p.multi     = [0,1,2]
  !p.background=truecolor('white')	; use a white background
  !p.color     =truecolor('black')	; use a black foreground
  
  ; plot the axis for emitted intensity as function of psi
  ttl = prefix+' - Beam emission intensity as function of ' + textoidl('\psi')
  plot, [0,1.0],[0,1.2*max(psi_emis)],/nodata,xs=1,ys=1,       $
        xtitle='normalised flux coordinate '+textoidl('\psi'), $
        ytitle='Intensity [photons/s/' + textoidl('\psi') + ']',         $
        title =ttl
  ; loop through the channels to plot the emitted intensity for each of them
  for k =0,nchan-1 do begin
    ; color the area underneath the emission white
    polyfill,[psi,reverse(psi)],$
             [psi_emis[*,k],replicate(0.0,n_elements(psi))],$
             color=truecolor('white')
    ; color the area that is used for the determination of the spatial resolution gray
    idx = where( (psi ge psi_res[3,k]) AND  (psi le psi_res[4,k]), count)
    if count ne 0 then begin
      polyfill,[psi[idx],reverse(psi[idx])], $
               [psi_emis[idx,k],replicate(0.0,count)], $
               color=truecolor('gray')
    endif
    ; plot the emission
    oplot, psi, psi_emis[*,k], thick=lth, linestyle=0
    ; plot the centre-of-mass psi
    mntmp = min( abs(psi-psi_res[2,k]), idx)
    oplot, [psi[idx],psi[idx]],[0,psi_emis[idx,k]],thick=1, linestyle=0
    ; plot the central psi
    mntmp = min( abs(psi-psi_res[0,k]), idx)
    oplot, [psi[idx],psi[idx]],[0,psi_emis[idx,k]],thick=1, linestyle=2
    ; plot the psi at maximum
    mntmp = min( abs(psi-psi_res[1,k]), idx)
    oplot, [psi[idx],psi[idx]],[0,psi_emis[idx,k]],thick=1, linestyle=4
  endfor
  ; plot the legend
  makelegend, ['centre of mass', 'central point', 'maximum'],$
          linestyle=[0,2,4], thick=[1,1,1], /box, /top, /left

  ; plot the axis for emitted intensity as function of R
  ttl = prefix+' - Beam emission intensity as function of R'
  plot, [min(R),max(R)],[0,1.2*max(R_emis)],/nodata,xs=1,ys=1,$
        xtitle='R [m]', ytitle='Intensity [photons/s/m]', title =ttl
  ; loop through the channels to plot the emitted intensity for each of them
  for k =0,nchan-1 do begin
    ; color the area underneath the emission white
    polyfill,[R,reverse(R)],$
             [R_emis[*,k],replicate(0.0,n_elements(R))],$
             color=truecolor('white')
    ; color the area that is used for the determination of the spatial resolution gray
    idx = where( (R ge R_res[3,k]) AND  (R le R_res[4,k]) , count)
    if count ne 0 then begin
      polyfill,[R[idx],reverse(R[idx])], $
               [R_emis[idx,k],replicate(0.0,n_elements(idx))], $
               color=truecolor('gray')
    endif
    ; plot the emission
    oplot, R, R_emis[*,k], thick=lth
    ; plot the centre-of-mass R
    mntmp = min( abs(R-R_res[2,k]), idx)
    oplot, [R[idx],R[idx]],[0,R_emis[idx,k]],thick=1, linestyle=0
    ; plot the central R
    mntmp = min( abs(R-R_res[0,k]), idx)
    oplot, [R[idx],R[idx]],[0,R_emis[idx,k]],thick=1, linestyle=2
    ; plot the R at maximum
    mntmp = min( abs(R-R_res[1,k]), idx)
    oplot, [R[idx],R[idx]],[0,R_emis[idx,k]],thick=1, linestyle=4
  endfor
  ; legend
  makelegend, ['centre of mass', 'central point', 'maximum'],$
          linestyle=[0,2,4], thick=[1,1,1], /box, /top, /left

  ; close the figure
  tmp = reset_graphic(figID, eps=epsfig)
  !p.multi = 0
endif


;----------------------------------------------------------------------------------
; plot the total emission intensity as function of psi and R if requested
;----------------------------------------------------------------------------------
if input.emisres.temisplot  then begin
  ; make the figures title
  ttl = prefix+' - Total emission intensity as function of '+ textoidl('\psi') + ' and R'

  ; set aspect ratio, figure zoom, fontsize and line thickness
  aspect   = 0.8
  figzoom  = 1.2
  fontsize = 12
  lth      = 2.0
 
  ; open a device for plotting the total emission as function of R and psi
  if (input.general.plot_print eq 1) then begin	
    figID = init_graphic(aspect=aspect, fontsize=fontsize, fontname=fontname, zoom=figzoom, title=ttl)
  endif
  if (input.general.plot_print eq 2) then begin	
    fname = pictdir+'/'+psprefix+'_total_emission_psi_R.eps'
    figID = init_graphic(aspect=aspect, fontsize=fontsize, fontname=fontname, /PSfont, zoom=figzoom, eps=fname, bits_per_pixel=8)
  endif
  !p.multi     = [0,1,2]
  !p.background=truecolor('white')	; use a white background
  !p.color     =truecolor('black')	; use a black foreground

  ; calculate the total emission intensity
  dR        = R[1]   - R[0]
  dpsi      = psi[1] - psi[0]
  R_temis   = total(R_emis,1)*dR
  psi_temis = total(psi_emis,1)*dpsi

  ; plot the total emission as function of psi
  ttl = prefix+' - Total emission intensity as function of ' + textoidl('\psi')
  xmin = 0
  xmax = 1
  ymin = 0
  ymax = 1.2*max(psi_temis)
  plot, psi_res[2,*], psi_temis, thick=lth, psym=-1,                                 $
        xs=1, xr=[xmin,xmax], ys=1, yr=[ymin,ymax],                                  $
        xtitle='normalised flux coordinate '+ textoidl('\psi') +' at Centre of Mass',$
        ytitle='intensity [photons/s]', title=ttl

  ; plot the total emission as function of R
  ttl = prefix+' - Total emission intensity as function of R'
  xmin = min(R)
  xmax = max(R)
  ymin = 0
  ymax = 1.2*max(R_temis)
  plot, R_res[2,*], R_temis, thick=lth, psym=-1,   $
        xs=1, xr=[xmin,xmax], ys=1, yr=[ymin,ymax],$
        xtitle='R - Centre of Mass [m]',           $
        ytitle='intensity [photons/s]', title=ttl

  ; close the figure
  tmp = reset_graphic(figID, eps=epsfig)
  !p.multi = 0
  ;----------------------------------------------------------------------------------
  ; print some information to the screen and the html-file
  ;----------------------------------------------------------------------------------
  print,''
  printf, funit,'<h2>Total emission intensity</h2>'
  printf, funit,'<p align="center"><table width="95%" border="0" cellspacing="0" cellpadding="0">'
  for k =0,nchan-1 do begin
    print,format='("  - total emitted intensity channel ",A,"  : ",(E9.2)," (from R), ",(E9.2)," (from psi)")',$
          chanID[k], total(R_emis[*,k])*dR, total(psi_emis[*,k])*dpsi
    printf,funit,'<tr><td width="37%" align="left">'
    printf,funit,format='("total emitted intensity channel ",A)',chanID[k]
    printf,funit,'</td><td width="5%" align="center">:</td><td width="58%" align="left">'
    printf,funit,format='((E9.2)," (from R), ",(E9.2)," (from psi)")',R_temis[k], psi_temis[k]
    printf,funit,'</td></tr>'
  endfor
  print,format='("  - average emitted intensity           : ",(E9.2)," (from R), ",(E9.2)," (from psi)")',$
        mean(R_temis), mean(psi_temis)

  printf,funit,'<tr height="7"><td width="37%" align="left"> </td><td width="5%" align="center"> </td>'
  printf,funit,'<td width="58%" align="left"> </td></tr>'
  printf,funit,'<tr><td width="37%" align="left">'
  printf,funit,format='("average emitted intensity")'
  printf,funit,'</td><td width="5%" align="center">:</td><td width="58%" align="left">'
  printf,funit,format='((E9.2)," (from R), ",(E9.2)," (from psi)")', mean(R_temis), mean(psi_temis)
  printf,funit,'</td></tr>'
  printf,funit,'</table></p>'
endif


;----------------------------------------------------------------------------------
; plot the spatial resolution as function of psi and R if requested
;----------------------------------------------------------------------------------
if input.emisres.resplot  then begin
  ; make the figures title
  ttl = prefix+' - Spatial resolution as function of '+ textoidl('\psi') + ' and R'

  ; set aspect ratio, figure zoom, fontsize and line thickness
  aspect   = 0.8
  figzoom  = 1.2
  fontsize = 12
  lth      = 2.0
 
  ; open a device for plotting the spatial resolution as function of R and psi
  if (input.general.plot_print eq 1) then begin	
    figID = init_graphic(aspect=aspect, fontsize=fontsize, fontname=fontname, zoom=figzoom, title=ttl)
  endif
  if (input.general.plot_print eq 2) then begin	
    fname = pictdir+'/'+psprefix+'_spat_resolution_psi_R.eps'
    figID = init_graphic(aspect=aspect, fontsize=fontsize, fontname=fontname, /PSfont, zoom=figzoom, eps=fname, bits_per_pixel=8)
  endif
  !p.multi     = [0,1,2]
  !p.background=truecolor('white')	; use a white background
  !p.color     =truecolor('black')	; use a black foreground

  ; plot the resolution as function of psi
  ttl = prefix+' - Spatial resolution as function of ' + textoidl('\psi')
  yttl= 'Spatial resolution ['+textoidl('\psi')+ '] - ' + string(format='(I2)',round(100*psi_res[6,0])) +'% emitted'
  xmin = 0
  xmax = 1
  ymin = 0
  ymax = 1.5*max(psi_res[5,*])
  plot, psi_res[2,*], psi_res[5,*], thick=lth, psym=-1,                              $
        xs=1, xr=[xmin,xmax], ys=1, yr=[ymin,ymax],                                  $
        xtitle='normalised flux coordinate '+ textoidl('\psi') +' at Centre of Mass',$
        ytitle=yttl, title=ttl

  ; plot the resolution as function of R
  ttl = prefix+' - Spatial resolution as function of R'
  yttl= 'Spatial resolution [m] - ' + string(format='(I2)',round(100*psi_res[6,0])) +'% emitted'
  xmin = min(R)
  xmax = max(R)
  ymin = 0
  ymax = 1.5*max(R_res[5,*])
  plot, R_res[2,*], R_res[5,*], thick=lth, psym=-1,$
        xs=1, xr=[xmin,xmax], ys=1, yr=[ymin,ymax],$
        xtitle='R - Centre of Mass [m]',           $
        ytitle=yttl, title=ttl

  ; close the figure
  tmp = reset_graphic(figID, eps=epsfig)
  !p.multi = 0

  ;----------------------------------------------------------------------------------
  ; print some information to the screen and the html-file
  ;----------------------------------------------------------------------------------
  print,''
  printf, funit,'<h2>Spatial resolution</h2>'
  printf, funit,'<p align="center"><table width="95%" border="0" cellspacing="0" cellpadding="0">'
  for k =0,nchan-1 do begin
    print,format='("  - spatial resolution channel ",A,"       : ",(F6.3)," (m)        , ",(F6.3)," (psi)")',$
          chanID[k], R_res[5,k], psi_res[5,k]

    printf,funit,'<tr><td width="37%" align="left">'
    printf,funit,format='("spatial resolution channel ",A)',chanID[k]
    printf,funit,'</td><td width="5%" align="center">:</td><td width="58%" align="left">'
    printf,funit,format='((F6.3)," (m)        , ",(F6.3)," (psi)")', R_res[5,k], psi_res[5,k]
    printf,funit,'</td></tr>'
  endfor
  print,format='("  - average spatial resolution          : ",(F6.3)," (m)        , ",(F6.3)," (psi)")',$
        total(R_res[5,*])/nchan, total(psi_res[5,*])/nchan

  printf,funit,'<tr height="7"><td width="37%" align="left"> </td><td width="5%" align="center"> </td>'
  printf,funit,'<td width="58%" align="left"> </td></tr>' ;<spacer height="10" type="block">
  printf,funit,'<tr><td width="37%" align="left">'
  printf,funit,format='("average spatial resolution")'
  printf,funit,'</td><td width="5%" align="center">:</td><td width="58%" align="left">'
  printf,funit,format='((F6.3)," (m)        , ",(F6.3)," (psi)")', total(R_res[5,*])/nchan, total(psi_res[5,*])/nchan
  printf,funit,'</td></tr>'
  printf,funit,'</table></p>'

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
    nplots=1.
    if (input.spec.multiplot) then begin
      nplots=0.
      if input.spec.spectrum then nplots++
      if input.spec.polangle then nplots++
      if input.spec.polfrac then nplots++

      ; make the figures title
      if R_res[2,k] lt 10 then ttl = string(format='(" - Channel ",A," at R=",(F4.2),"m / '+textoidl('\psi')+'=",(F4.2))',$
                                     chanID[k],R_res[2,k],psi_res[2,k])                                                   $
                          else ttl = string(format='(" - Channel ",A," at R=",(I0),"m")',                                 $
                                     chanID[k],R_res[2,k])
      ttl = prefix+ttl

      ; set aspect ratio, figure zoom, fontsize and line thickness
      aspect   = 2.5/sqrt(nplots*3.)
      figzoom  = 1.2
      fontsize = 12*sqrt(nplots)
      lth      = 2.0
	 
      ; open a device for plotting the spectra as function of R and psi
      if (input.general.plot_print eq 1) then begin	
        figID = init_graphic(aspect=aspect, fontsize=fontsize, fontname=fontname, zoom=figzoom, title=ttl)
      endif
      if (input.general.plot_print eq 2) then begin	
        if R_res[2,k] lt 10 then fname = string(format='("_chan",A,"_R",(F4.2),"_psi",(F4.2),".eps")',$
                                                chanID[k],R_res[2,k],psi_res[2,k])                    $
                            else fname = string(format='("_chan",A,"_R",(I0),".eps")',                $
                                                chanID[k],R_res[2,k])
        fname = pictdir+'/'+psprefix+fname
        figID = init_graphic(aspect=aspect, fontsize=fontsize, fontname=fontname, /PSfont, zoom=figzoom, eps=fname, bits_per_pixel=8)
      endif
      !p.multi=[0,1,nplots]
      !p.background=truecolor('white')	; use a white background
      !p.color     =truecolor('black')	; use a black foreground
      !x.margin    = [10,9]
    endif


    ;----------------------------------------------------------------------------------
    ; Convert the Stokes vectors in (nonzero) intensity, polarised fraction,
    ;  polarisation angle and S/N 
    ;----------------------------------------------------------------------------------
    nonzero  = where(stokes[*,0,k] gt  0.0,count)   ; where the total intensity isn't zero
    tspec  = fltarr(nlambda)
    tpspec = fltarr(nlambda)
    ppspec = fltarr(nlambda)
    puspec = fltarr(nlambda)
    spspec = fltarr(nlambda)
    suspec = fltarr(nlambda)
    tpfrac = fltarr(nlambda)+1.0
    SN     = fltarr(nlambda)
    pgamma = fltarr(nlambda)
    sgamma = fltarr(nlambda)
    gamma  = fltarr(nlambda)
    if count ne 0 then begin
      tspec[nonzero]  = stokes[nonzero,0,k]
      tpspec[nonzero] = sqrt(stokes[nonzero,1,k]^2+stokes[nonzero,2,k]^2) ;Total linear polarization intensity

      ppspec[nonzero] = sqrt(pstokes[nonzero,1,k]^2+pstokes[nonzero,2,k]^2) ;Linear polarization intensity of the pi-lines
      puspec[nonzero] = pstokes[nonzero,0,k] - sqrt(pstokes[nonzero,1,k]^2+pstokes[nonzero,2,k]^2) ; Circular polarization intensity + unpolarized intensity of the pi-lines

      spspec[nonzero] = sqrt(sstokes[nonzero,1,k]^2+sstokes[nonzero,2,k]^2) ;sigma linear polarization intensity
      suspec[nonzero] = sstokes[nonzero,0,k] - sqrt(sstokes[nonzero,1,k]^2+sstokes[nonzero,2,k]^2) ;sigma unpolarised
      tpfrac[nonzero] = tpspec[nonzero]/tspec[nonzero] ;total polarised fraction
      SN[nonzero]     = tpspec[nonzero]/sqrt(tspec[nonzero]) ;signal to noise
      
      tpnonzero = where(tpspec gt 0.0)  ; where the total polarised intensity,
      ppnonzero = where(ppspec gt 0.0)  ; the polarised pi-intensity
      spnonzero = where(spspec gt 0.0)  ; and the polarised sigma-intensity isn't zero

      pgamma[ppnonzero] = 0.5*atan(pstokes[ppnonzero,2,k],pstokes[ppnonzero,1,k])*!radeg
      sgamma[spnonzero] = 0.5*atan(sstokes[spnonzero,2,k],sstokes[spnonzero,1,k])*!radeg
      gamma[tpnonzero]  = 0.5*atan(stokes[tpnonzero,2,k],stokes[tpnonzero,1,k])*!radeg
    endif

    ;----------------------------------------------------------------------------------
    ; Set the wavelength limits for the plotting
    ;----------------------------------------------------------------------------------
    range  = where(stokes[*,0,k]-min(stokes[*,0,k]) gt 0.05*(max(stokes[*,0,k])-min(stokes[*,0,k])),nrange) ; where the intensity is larger than 
    if nrange ne 0 then begin                                                                               ; 5% of the difference between maximum
      lambdamin = lambda[range[0]]-4.0                                                                      ; and minimum intensity
      lambdamax = lambda[range[nrange-1]]+4.0
    endif else begin
      lambdamin = lambda0+Dshift0[k]-4.0
      lambdamax = lambda0+Dshift0[k]+4.0
    endelse


    ;----------------------------------------------------------------------------------
    ; Locate the CWL and the wavelengths of max. pi-blue, max. sigma and max.pi-red
    ;----------------------------------------------------------------------------------
    tmp = min(abs(lambda-cwlstokes[4,k]),cwlidx)
    tmp = min(abs(lambda-pibstokes[4,k]),pbidx)
    tmp = min(abs(lambda-sigstokes[4,k]),sidx)
    tmp = min(abs(lambda-pirstokes[4,k]),pridx)
    
    ;----------------------------------------------------------------------------------
    ; print some information to the screen and the html-file
    ;----------------------------------------------------------------------------------
    print,''
    print,format='("  - Spectral information for channel ",A," at R=",F4.2,"m")',chanID[k],R_res[2,k]
    printf, funit,format='("<h2>Spectral information for channel ",A," at R=",F4.2,"m</h2>")',chanID[k],R_res[2,k]

    ;----------------------------------------------------------------------------------
    ; plot spectrum if requested 
    ;----------------------------------------------------------------------------------
    if input.spec.spectrum then begin

      ;----------------------------------------------------------------------------------
      ; If multiplot is NOT selected, then we open a device for plotting here.
      ;----------------------------------------------------------------------------------
      if ~(input.spec.multiplot) then begin
        ; set aspect ratio, figure zoom, fontsize and line thickness
        aspect   = 1.5
        figzoom  = 1.2
        fontsize = 12
        lth      = 2.0
		 
        ; open a device for plotting the spectra as function of wavelength
        if (input.general.plot_print eq 1) then begin	
        figID = init_graphic(aspect=aspect, fontsize=fontsize, fontname=fontname, zoom=figzoom)
        endif
        if (input.general.plot_print eq 2) then begin	
        if R_res[2,k] lt 10 then fname = string(format='("_spec",A,"_R",(F4.2),"_psi",(F4.2),".eps")',$
                                                chanID[k],R_res[2,k],psi_res[2,k])                    $
                            else fname = string(format='("_spec",A,"_R",(I0),".eps")',                $
                                                chanID[k],R_res[2,k])
          fname = pictdir+'/'+psprefix+fname
          figID = init_graphic(aspect=aspect, fontsize=fontsize, fontname=fontname, /PSfont, zoom=figzoom, eps=fname, bits_per_pixel=8)
        endif
        !p.multi     = 0
        !p.background= truecolor('white')	; use a white background
        !p.color     = truecolor('black')	; use a black foreground
      endif

      ;----------------------------------------------------------------------------------
      ; make the figures title
      ;----------------------------------------------------------------------------------
      if R_res[2,k] lt 10 then ttl = string(format='(" - Spectrum ",A," at R=",(F4.2),"m / '+textoidl('\psi')+'=",(F4.2))',$
                                     chanID[k],R_res[2,k],psi_res[2,k])                                                    $
                          else ttl = string(format='(" - Spectrum ",A," at R=",(I0),"m")',                                 $
                                     chanID[k],R_res[2,k])
      ttl = prefix+ttl

      ;----------------------------------------------------------------------------------
      ; the plotting
      ;----------------------------------------------------------------------------------
      ; set y-limits
      ymin = 0.
      ymax = 1.05 * max(tspec)
      ; total spectrum
      if filterflag eq 0 then yttl='Intensity [photons/s/Angstrom]' else yttl='Intensity [photons/s]'
      if filterflag eq 0 then xttl='Wavelength [Angstrom]' else xttl='Filter central wavelength [Angstrom]'
      plot, lambda, tspec, color=truecolor('black'), thick=lth, linestyle=0,$
            xs=1, xr=[lambdamin,lambdamax], ys=1, yr=[ymin,ymax],           $
            xtitle=xttl, ytitle=yttl, title=ttl
      ; pol. pi spectrum
      oplot,lambda, ppspec, color=truecolor('blue'), thick=lth, linestyle=2
      ; pol. sigma spectrum
      oplot,lambda, spspec, color=truecolor('red'), thick=lth, linestyle=3
      ; unpol. pi spectrum
      oplot,lambda, puspec, color=truecolor('darkblue'), thick=lth, linestyle=4
      ; unpol. sigma spectrum
      oplot,lambda, suspec, color=truecolor('darkred'), thick=lth, linestyle=5
      ; indicate the CWL, max. S/N blue shifted pi, red shifted pi and sigma
      oplot, [cwlstokes[4,k],cwlstokes[4,k]],[ymin,ymax],color=truecolor('black'),thick=1,linestyle=1
      plots, cwlstokes[4,k], cwlstokes[0,k],psym=1,color=truecolor('black'),thick=1, /data
      oplot, [pibstokes[4,k],pibstokes[4,k]],[ymin,ymax],color=truecolor('blue'),thick=1,linestyle=1
      plots, pibstokes[4,k], pibstokes[0,k],psym=1,color=truecolor('blue'),thick=1, /data
      oplot, [pirstokes[4,k],pirstokes[4,k]],[ymin,ymax],color=truecolor('blue'),thick=1,linestyle=1
      plots, pirstokes[4,k], pirstokes[0,k],psym=1,color=truecolor('blue'),thick=1, /data
      oplot, [sigstokes[4,k],sigstokes[4,k]],[ymin,ymax],color=truecolor('red'),thick=1,linestyle=1
      plots, sigstokes[4,k], sigstokes[0,k],psym=1,color=truecolor('red'),thick=1, /data
      ; legend
      makelegend, ['total intensity',                                                                     $
               'linear pol.' + textoidl('\pi'),    'circular and upol.' + textoidl('\pi'),                                        $
               'linear pol.' + textoidl('\sigma'), 'circular and upol.' + textoidl('\sigma'),                                     $
               textoidl('\lambda') +'cwl',                                                             $
               textoidl('\lambda') +'max pi',                                                         $
               textoidl('\lambda') + 'max sigma'],                                                    $
               textcolors=truecolor(['black','blue','darkblue','red','darkred','black','blue','red']),$
               colors    =truecolor(['black','blue','darkblue','red','darkred','black','blue','red']),$
               linestyle=[0,2,3,4,5,1,1,1],thick=[lth,lth,lth,lth,lth,1,1,1],                         $
               charsize=1./sqrt(nplots), spacing=1./sqrt(nplots), /clear, /box, /top, /left

      ;----------------------------------------------------------------------------------
      ; close the device if needed
      ;----------------------------------------------------------------------------------
      if ~(input.spec.multiplot) then begin
        ; close the figure
        tmp = reset_graphic(figID, eps=epsfig)
        !p.multi = 0
      endif

      ;----------------------------------------------------------------------------------
      ; print some information to the screen and the html-file
      ;----------------------------------------------------------------------------------
      if filterflag eq 0 then unit='photons/s/A' else unit='photons/s'
      if filterflag eq 0 then cwltxt='of the spectrum' else cwltxt='of the filter  '

      print,format='("     - central wavelength (CWL) ",A,": ", (F7.2)," A")',$
            cwltxt, lambda[cwlidx]
      print,format='("     - wavelength at max. pi-blue                : ", (F7.2)," A")',$
            lambda[pbidx]
      print,format='("     - wavelength at max. sigma                  : ", (F7.2)," A")',$
            lambda[sidx]
      print,format='("     - wavelength at max. pi-red                 : ", (F7.2)," A")',$
            lambda[pridx]
      print,''

      printf,funit,'<p align="center"><table width="95%" border="0" cellspacing="0" cellpadding="0">'
      printf,funit,'<tr><td width="37%" align="left">'
      printf,funit,format='("central wavelength (CWL) ",A)', cwltxt
      printf,funit,'</td>'
      printf,funit,'<td width="5%" align="center">:</td><td width="58%" align="left">'
      printf,funit,format='(F7.2," A")', lambda[cwlidx]
      printf,funit,'</td></tr>'
      printf,funit,'<tr><td width="37%" align="left">wavelength at max. pi-blue</td>'
      printf,funit,'<td width="5%" align="center">:</td><td width="58%" align="left">'
      printf,funit,format='(F7.2," A")', lambda[pbidx]
      printf,funit,'</td></tr>'
      printf,funit,'<tr><td width="37%" align="left">wavelength at max. sigma</td>'
      printf,funit,'<td width="5%" align="center">:</td><td width="58%" align="left">'
      printf,funit,format='(F7.2," A")', lambda[sidx]
      printf,funit,'</td></tr>'
      printf,funit,'<tr><td width="37%" align="left">wavelength at max. pi-red</td>'
      printf,funit,'<td width="5%" align="center">:</td><td width="58%" align="left">'
      printf,funit,format='(F7.2," A")', lambda[pridx]
      printf,funit,'</td></tr>'
      printf,funit,'</table></p>'

      print,format='("     - total intensity at the CWL                : ",(E10.1)," ",A)',$
            tspec[cwlidx], unit
      print,format='("     - total intensity at max. pi-blue           : ",(E10.1)," ",A)',$
            tspec[pbidx], unit
      print,format='("     - total intensity at max. sigma             : ",(E10.1)," ",A)',$
            tspec[sidx], unit
      print,format='("     - total intensity at max. pi-red            : ",(E10.1)," ",A)',$
            tspec[pridx], unit
      print,''

      printf,funit,'<p align="center"><table width="95%" border="0" cellspacing="0" cellpadding="0">'
      printf,funit,'<tr><td width="37%" align="left">total intensity at the CWL</td>'
      printf,funit,'<td width="5%" align="center">:</td><td width="58%" align="left">'
      printf,funit,format='((E10.1)," ",A)', tspec[cwlidx], unit
      printf,funit,'</td></tr>'
      printf,funit,'<tr><td width="37%" align="left">total intensity at max. pi-blue</td>'
      printf,funit,'<td width="5%" align="center">:</td><td width="58%" align="left">'
      printf,funit,format='((E10.1)," ",A)', tspec[pbidx], unit
      printf,funit,'</td></tr>'
      printf,funit,'<tr><td width="37%" align="left">total intensity at max. sigma</td>'
      printf,funit,'<td width="5%" align="center">:</td><td width="58%" align="left">'
      printf,funit,format='((E10.1)," ",A)', tspec[sidx], unit
      printf,funit,'</td></tr>'
      printf,funit,'<tr><td width="37%" align="left">total intensity at max. pi-red</td>'
      printf,funit,'<td width="5%" align="center">:</td><td width="58%" align="left">'
      printf,funit,format='((E10.1)," ",A)', tspec[pridx], unit
      printf,funit,'</td></tr>'
      printf,funit,'</table></p>'

    endif


    ;----------------------------------------------------------------------------------
    ; plot pol. angle if requested 
    ;----------------------------------------------------------------------------------
    if input.spec.polangle then begin

      ;----------------------------------------------------------------------------------
      ; If multiplot is NOT selected, then we open a device for plotting here.
      ;----------------------------------------------------------------------------------
      if ~(input.spec.multiplot) then begin
        ; set aspect ratio, figure zoom, fontsize and line thickness
        aspect   = 1.5
        figzoom  = 1.2
        fontsize = 12
        lth      = 2.0
		 
        ; open a device for plotting the spectra as function of wavelength
        if (input.general.plot_print eq 1) then begin	
          figID = init_graphic(aspect=aspect, fontsize=fontsize, fontname=fontname, zoom=figzoom)
        endif
        if (input.general.plot_print eq 2) then begin	
        if R_res[2,k] lt 10 then fname = string(format='("_polangle",A,"_R",(F4.2),"_psi",(F4.2),".eps")',$
                                                chanID[k],R_res[2,k],psi_res[2,k])                    $
                            else fname = string(format='("_polangle",A,"_R",(I0),".eps")',                $
                                                chanID[k],R_res[2,k])
          fname = pictdir+'/'+psprefix+fname
          figID = init_graphic(aspect=aspect, fontsize=fontsize, fontname=fontname, /PSfont, zoom=figzoom, eps=fname, bits_per_pixel=8)
        endif
        !p.multi     = 0
        !p.background= truecolor('white')	; use a white background
        !p.color     = truecolor('black')	; use a black foreground
      endif

      ;----------------------------------------------------------------------------------
      ; make the figures title
      ;----------------------------------------------------------------------------------
      if R_res[2,k] lt 10 then ttl = string(format='(" - Polarisation Angle ",A," at R=",(F4.2),"m / '+textoidl('\psi')+'=",(F4.2))',$
                                     chanID[k],R_res[2,k],psi_res[2,k])                                                              $
                          else ttl = string(format='(" - Polarisation Angle ",A," at R=",(I0),"m")',                                 $
                                     chanID[k],R_res[2,k])
      ttl = prefix+ttl

      ;----------------------------------------------------------------------------------
      ; the plotting of the pol. angles (only where the polarised intensity is nonzero)
      ;----------------------------------------------------------------------------------
      ; set y-limits
      ymin = -90.0
      ymax = 90.0
      ; overall polarisation angle
      if filterflag eq 0 then xttl='Wavelength [Angstrom]' else xttl='Filter central wavelength [Angstrom]'
      plot, lambda[tpnonzero], gamma[tpnonzero], thick=lth, linestyle=0,$
            xs=1, xr=[lambdamin,lambdamax], ys=1, yr=[ymin,ymax],$
            xtitle=xttl, ytitle='Pol. angle '+textoidl('\gamma')+' [degrees]', $
            title=ttl
      ; pi-polarisation angle
      oplot, lambda[ppnonzero], pgamma[ppnonzero], color=truecolor('blue'), thick=lth, linestyle=2
      ; sigma-polarisation angle
      oplot, lambda[spnonzero], sgamma[spnonzero], color=truecolor('red'), thick=lth, linestyle=2
      ; plot ideal polarization angle at 0, -90, -180
      oplot, [lambdamin,lambdamax], [alpha0[k]*!radeg,alpha0[k]*!radeg], color=truecolor('black'), thick=0.5, linestyle=3 
      oplot, [lambdamin,lambdamax], [alpha0[k]*!radeg-90,alpha0[k]*!radeg-90], color=truecolor('black'), thick=0.5, linestyle=3 
      oplot, [lambdamin,lambdamax], [alpha0[k]*!radeg-180,alpha0[k]*!radeg-180], color=truecolor('black'), thick=0.5, linestyle=3 
      ; indicate the CWL, max. S/N blue shifted pi, red shifted pi and sigma
      oplot, [cwlstokes[4,k],cwlstokes[4,k]],[ymin,ymax],color=truecolor('black'),thick=1,linestyle=1
      plots, cwlstokes[4,k], 0.5*atan(cwlstokes[2,k],cwlstokes[1,k])*!radeg,psym=1,color=truecolor('black'),thick=1, /data
      oplot, [pibstokes[4,k],pibstokes[4,k]],[ymin,ymax],color=truecolor('blue'),thick=1,linestyle=1
      plots, pibstokes[4,k], 0.5*atan(pibstokes[2,k],pibstokes[1,k])*!radeg,psym=1,color=truecolor('blue'),thick=1, /data
      oplot, [pirstokes[4,k],pirstokes[4,k]],[ymin,ymax],color=truecolor('blue'),thick=1,linestyle=1
      plots, pirstokes[4,k], 0.5*atan(pirstokes[2,k],pirstokes[1,k])*!radeg,psym=1,color=truecolor('blue'),thick=1, /data
      oplot, [sigstokes[4,k],sigstokes[4,k]],[ymin,ymax],color=truecolor('red'),thick=1,linestyle=1
      plots, sigstokes[4,k], 0.5*atan(sigstokes[2,k],sigstokes[1,k])*!radeg,psym=1,color=truecolor('red'),thick=1, /data
      ; legend
      makelegend, ['total polarisation angle' + textoidl('\gamma'),                                  $
               textoidl('\pi') + 'polarisation angle' + textoidl('\gamma \pi'),                                $
               textoidl('\sigma') + 'polarisation angle' + textoidl('\gamma_\sigma'),                          $
               'ideal polarisation angle' + textoidl('\gamma'),                        $
               textoidl('\lambda') +'cwl',                                                    $
               textoidl('\lambda') +'max pi',                                                $
               textoidl('\lambda') + 'max sigma'],                                           $
               textcolors=truecolor(['black','blue','red','black','black','blue','red']),    $
               colors    =truecolor(['black','blue','red','black','black','blue','red']),    $
               linestyle=[0,2,2,3,1,1,1],thick=[lth,lth,lth,0.5,1,1,1],                      $
               charsize=1./sqrt(nplots), spacing=1./sqrt(nplots), /clear, /box, /top, /left

      ;----------------------------------------------------------------------------------
      ; close the device if needed
      ;----------------------------------------------------------------------------------
      if ~(input.spec.multiplot) then begin
        ; close the figure
        tmp = reset_graphic(figID, eps=epsfig)
        !p.multi = 0
      endif

      ;----------------------------------------------------------------------------------
      ; print some information to the screen and the html-file
      ;----------------------------------------------------------------------------------
      cwlrange = [(cwlidx-5)>0,(cwlidx+5)<(nlambda-1)] ; range over which the change in polarisation angle
      pbrange  = [(pbidx-5)>0 ,(pbidx+5)<(nlambda-1) ] ; and the change in S/N-ratio is determined
      srange   = [(sidx-5)>0  ,(sidx+5)<(nlambda-1)  ]
      prrange  = [(pridx-5)>0 ,(pridx+5)<(nlambda-1) ]
      print,format='("     - pol. angle at the CWL                     : ",(F6.2)," degrees")',$
            gamma[cwlidx]
      print,format='("       change in pol. angle at the CWL           : ",(F6.2)," degrees/Angstrom")',$
            (gamma[cwlrange[1]] - gamma[cwlrange[0]])/(lambda[cwlrange[1]] - lambda[cwlrange[0]])
      print,format='("     - pol. angle at max. pi-blue                : ",(F6.2)," degrees")',$
            gamma[pbidx]
      print,format='("       change in pol. angle at pi-blue           : ",(F6.2)," degrees/Angstrom")',$
            (gamma[pbrange[1]] - gamma[pbrange[0]])/(lambda[pbrange[1]] - lambda[pbrange[0]])
      print,format='("     - pol. angle at max. sigma                  : ",(F6.2)," degrees")',$
            gamma[sidx]
      print,format='("       change in pol. angle at sigma             : ",(F6.2)," degrees/Angstrom")',$
            (gamma[srange[1]] - gamma[srange[0]])/(lambda[srange[1]] - lambda[srange[0]])
      print,format='("     - pol. angle at max. pi-red                 : ",(F6.2)," degrees")',$
            gamma[pridx]
      print,format='("       change in pol. angle at pi-red            : ",(F6.2)," degrees/Angstrom")',$
            (gamma[prrange[1]] - gamma[prrange[0]])/(lambda[prrange[1]] - lambda[prrange[0]])
      print,''

      printf, funit,'<p align="center"><table width="95%" border="0" cellspacing="0" cellpadding="0">'
      printf,funit,'<tr><td width="37%" align="left">pol. angle at the CWL</td>'
      printf,funit,'<td width="5%" align="center">:</td><td width="58%" align="left">'
      printf,funit,format='((F6.2)," degrees")', gamma[cwlidx]
      printf,funit,'</td></tr>'
      printf,funit,'<tr><td width="37%" align="left">change in pol. angle at the CWL</td>'
      printf,funit,'<td width="5%" align="center">:</td><td width="58%" align="left">'
      printf,funit,format='((F6.2)," degrees/Angstrom")',(gamma[cwlrange[1]]-gamma[cwlrange[0]])/(lambda[cwlrange[1]]-lambda[cwlrange[0]])
      printf,funit,'</td></tr>'
      printf,funit,'<tr><td width="37%" align="left">pol. angle at max. pi-blue</td>'
      printf,funit,'<td width="5%" align="center">:</td><td width="58%" align="left">'
      printf,funit,format='((F6.2)," degrees")', gamma[pbidx]
      printf,funit,'</td></tr>'
      printf,funit,'<tr><td width="37%" align="left">change in pol. angle at pi-blue</td>'
      printf,funit,'<td width="5%" align="center">:</td><td width="58%" align="left">'
      printf,funit,format='((F6.2)," degrees/Angstrom")',(gamma[pbrange[1]]-gamma[pbrange[0]])/(lambda[pbrange[1]]-lambda[pbrange[0]])
      printf,funit,'</td></tr>'
      printf,funit,'<tr><td width="37%" align="left">pol. angle at max. sigma</td>'
      printf,funit,'<td width="5%" align="center">:</td><td width="58%" align="left">'
      printf,funit,format='((F6.2)," degrees")', gamma[sidx]
      printf,funit,'</td></tr>'
      printf,funit,'<tr><td width="37%" align="left">change in pol. angle at sigma</td>'
      printf,funit,'<td width="5%" align="center">:</td><td width="58%" align="left">'
      printf,funit,format='((F6.2)," degrees/Angstrom")',(gamma[srange[1]]-gamma[srange[0]])/(lambda[srange[1]]-lambda[srange[0]])
      printf,funit,'</td></tr>'
      printf,funit,'<tr><td width="37%" align="left">pol. angle at max. pi-red</td>'
      printf,funit,'<td width="5%" align="center">:</td><td width="58%" align="left">'
      printf,funit,format='((F6.2)," degrees")', gamma[pridx]
      printf,funit,'</td></tr>'
      printf,funit,'<tr><td width="37%" align="left">change in pol. angle at pi-red</td>'
      printf,funit,'<td width="5%" align="center">:</td><td width="58%" align="left">'
      printf,funit,format='((F6.2)," degrees/Angstrom")',(gamma[prrange[1]]-gamma[prrange[0]])/(lambda[prrange[1]]-lambda[prrange[0]])
      printf,funit,'</td></tr>'
      printf,funit,'</table></p>'

    endif

    ;----------------------------------------------------------------------------------
    ; plot pol. fraction and norm. S/N ratio if requested 
    ;----------------------------------------------------------------------------------
    if input.spec.polfrac then begin

      ;----------------------------------------------------------------------------------
      ; If multiplot is NOT selected, then we open a device for plotting here.
      ;----------------------------------------------------------------------------------
      if ~(input.spec.multiplot) then begin
        ; set aspect ratio, figure zoom, fontsize and line thickness
        aspect   = 1.5
        figzoom  = 1.2
        fontsize = 12
        lth      = 2.0
		 
        ; open a device for plotting the spectra as function of wavelength
        if (input.general.plot_print eq 1) then begin	
          figID = init_graphic(aspect=aspect, fontsize=fontsize, fontname=fontname, zoom=figzoom)
        endif
        if (input.general.plot_print eq 2) then begin	
        if R_res[2,k] lt 10 then fname = string(format='("_polfrac",A,"_R",(F4.2),"_psi",(F4.2),".eps")',$
                                                chanID[k],R_res[2,k],psi_res[2,k])                    $
                            else fname = string(format='("_polfrac",A,"_R",(I0),".eps")',                $
                                                chanID[k],R_res[2,k])
          fname = pictdir+'/'+psprefix+fname
          figID = init_graphic(aspect=aspect, fontsize=fontsize, fontname=fontname, /PSfont, zoom=figzoom, eps=fname, bits_per_pixel=8)
        endif
        !p.multi     = 0
        !p.background= truecolor('white')	; use a white background
        !p.color     = truecolor('black')	; use a black foreground
        !x.margin    = [10,9]
      endif

      ;----------------------------------------------------------------------------------
      ; make the figures title
      ;----------------------------------------------------------------------------------
      if R_res[2,k] lt 10 then ttl = string(format='(" - Pol. fraction/SN-ratio ",A," at R=",(F4.2),"m / '+textoidl('\psi')+'=",(F4.2))',$
                                     chanID[k],R_res[2,k],psi_res[2,k])                                                                  $
                          else ttl = string(format='(" - Pol. fraction/SN-ratio ",A," at R=",(I0),"m")',                                 $
                                     chanID[k],R_res[2,k])
      ttl = prefix+ttl

      ;----------------------------------------------------------------------------------
      ; the plotting of the pol. fraction (only where the polarised intensity is nonzero)
      ;----------------------------------------------------------------------------------
      ; set y-limits
      ymin = 0.0
      ymax = 1.0
      ; polarised fraction
      if filterflag eq 0 then xttl='Wavelength [Angstrom]' else xttl='Filter central wavelength [Angstrom]'
      plot, lambda[tpnonzero], tpfrac[tpnonzero], thick=lth,$
            linestyle=0, xs=1, xr=[lambdamin,lambdamax], ys=5, yr=[ymin,ymax],$
            xtitle=xttl, $
            title=ttl
      axis, yaxis=0, ys=1, yr=[ymin,ymax], ytitle='Linear polarisation fraction'
      ; S/N ratio
      position = [!x.window, !y.window]
      position = position[[0,2,1,3]]
      ymin=0
      ymax=1.05*max(SN)
      if filterflag eq 0 then yttl='S/N-ratio [(photons/s/Angstrom)^{1/2}]' $
                         else yttl='S/N-ratio [(photons/s)^{1/2}]'
      plot, lambda[tpnonzero], SN[tpnonzero], linestyle=2, thick=lth, $
            xs=5, xr=[lambdamin,lambdamax], ys=5, yr=[ymin,ymax],     $
            position=position, /noerase
      axis, yaxis=1, ys=1, yr=[ymin,ymax], ytitle=yttl
      ; indicate the CWL, max. S/N blue shifted pi, red shifted pi and sigma
      oplot, [cwlstokes[4,k],cwlstokes[4,k]],[ymin,ymax],color=truecolor('black'),thick=1,linestyle=1
      oplot, [pibstokes[4,k],pibstokes[4,k]],[ymin,ymax],color=truecolor('blue'),thick=1,linestyle=1
      oplot, [pirstokes[4,k],pirstokes[4,k]],[ymin,ymax],color=truecolor('blue'),thick=1,linestyle=1
      oplot, [sigstokes[4,k],sigstokes[4,k]],[ymin,ymax],color=truecolor('red'),thick=1,linestyle=1
      ; legend
      makelegend, ['linear polarisation fraction',                                     $
               'S/N ratio',                                                        $
               textoidl('\lambda')+'cwl',                                                    $
               textoidl('\lambda')+'max pi',                                                $
               textoidl('\lambda')+ 'max sigma'],                                           $
               textcolors=truecolor(['black','black','black','blue','red']),                 $
               colors    =truecolor(['black','black','black','blue','red']),                 $
               linestyle=[0,2,1,1,1],thick=[lth,lth,1,1,1],                                  $
               charsize=1./sqrt(nplots), spacing=1./sqrt(nplots), /clear, /box, /top, /left

      ;----------------------------------------------------------------------------------
      ; close the device if needed
      ;----------------------------------------------------------------------------------
      if ~(input.spec.multiplot) then begin
        ; close the figure
        tmp = reset_graphic(figID, eps=epsfig)
        !p.multi = 0
      endif

      ;----------------------------------------------------------------------------------
      ; print some information to the screen and the html-file
      ;----------------------------------------------------------------------------------
      if filterflag eq 0 then unit='(photons/s/A)^{1/2}' else unit='(photons/s)^{1/2}'
      print,format='("     - pol. fraction at the CWL                  : ",(I6)," %")',$
            round(100*tpfrac[cwlidx])
      print,format='("     - pol. fraction at max. pi-blue             : ",(I6)," %")',$
            round(100*tpfrac[pbidx])
      print,format='("     - pol. fraction at max. sigma               : ",(I6)," %")',$
            round(100*tpfrac[sidx])
      print,format='("     - pol. fraction at max. pi-red              : ",(I6)," %")',$
            round(100*tpfrac[pridx])
      print,''

      printf, funit,'<p align="center"><table width="95%" border="0" cellspacing="0" cellpadding="0">'
      printf,funit,'<tr><td width="37%" align="left">pol. fraction at the cwl</td>'
      printf,funit,'<td width="5%" align="center">:</td><td width="58%" align="left">'
      printf,funit,format='((I6)," %")', round(100*tpfrac[cwlidx])
      printf,funit,'</td></tr>'
      printf,funit,'<tr><td width="37%" align="left">pol. fraction at max. pi-blue</td>'
      printf,funit,'<td width="5%" align="center">:</td><td width="58%" align="left">'
      printf,funit,format='((I6)," %")', round(100*tpfrac[pbidx])
      printf,funit,'</td></tr>'
      printf,funit,'<tr><td width="37%" align="left">pol. fraction at max. sigma</td>'
      printf,funit,'<td width="5%" align="center">:</td><td width="58%" align="left">'
      printf,funit,format='((I6)," %")', round(100*tpfrac[sidx])
      printf,funit,'</td></tr>'
      printf,funit,'<tr><td width="37%" align="left">pol. fraction at max. pi-red</td>'
      printf,funit,'<td width="5%" align="center">:</td><td width="58%" align="left">'
      printf,funit,format='((I6)," %")', round(100*tpfrac[pridx])
      printf,funit,'</td></tr>'
      printf,funit,'</table></p>'
    
      print,format='("     - S/N ratio at the CWL                      : ",E10.2," ", A," (",I3," % of maximum)")',$
            SN[cwlidx], unit, round(100*SN[cwlidx]/max(SN))
      print,format='("       CWL S/N ratio at +",F4.2,"A                   : ",E10.2," ", A," (",I3," % of maximum)")',$
            lambda[cwlrange[1]]-lambda[cwlidx], SN[cwlrange[1]], unit, round(100*SN[cwlrange[1]]/max(SN))
      print,format='("       CWL S/N ratio at -",F4.2,"A                   : ",E10.2," ", A," (",I3," % of maximum)")',$
            lambda[cwlidx]-lambda[cwlrange[0]], SN[cwlrange[0]], unit, round(100*SN[cwlrange[0]]/max(SN))
      print,''
      print,format='("     - S/N ratio at max. pi-blue                 : ",E10.2," ", A," (",I3," % of maximum)")',$
            SN[pbidx], unit, round(100*SN[pbidx]/max(SN))
      print,format='("       pi-blue S/N ratio at +",F4.2,"A               : ",E10.2," ", A," (",I3," % of maximum)")',$
            lambda[pbrange[1]]-lambda[pbidx], SN[pbrange[1]], unit, round(100*SN[pbrange[1]]/max(SN))
      print,format='("       pi-blue S/N ratio at -",F4.2,"A               : ",E10.2," ", A," (",I3," % of maximum)")',$
            lambda[pbidx]-lambda[pbrange[0]], SN[pbrange[0]], unit, round(100*SN[pbrange[0]]/max(SN))
      print,''
      print,format='("     - S/N ratio at max. sigma                   : ",E10.2," ", A," (",I3," % of maximum)")',$
            SN[sidx], unit, round(100*SN[sidx]/max(SN))
      print,format='("       sigma S/N ratio at +",F4.2,"A                 : ",E10.2," ", A," (",I3," % of maximum)")',$
            lambda[srange[1]]-lambda[sidx], SN[srange[1]], unit, round(100*SN[srange[1]]/max(SN))
      print,format='("       sigma S/N ratio at -",F4.2,"A                 : ",E10.2," ", A," (",I3," % of maximum)")',$
            lambda[sidx]-lambda[srange[0]], SN[srange[0]], unit, round(100*SN[srange[0]]/max(SN))
      print,''
      print,format='("     - S/N ratio at max. pi-red                  : ",E10.2," ", A," (",I3," % of maximum)")',$
            SN[pridx], unit, round(100*SN[pridx]/max(SN))
      print,format='("       pi-red S/N ratio at +",F4.2,"A                : ",E10.2," ", A," (",I3," % of maximum)")',$
            lambda[prrange[1]]-lambda[pridx], SN[prrange[1]], unit, round(100*SN[prrange[1]]/max(SN))
      print,format='("       pi-red S/N ratio at -",F4.2,"A                : ",E10.2," ", A," (",I3," % of maximum)")',$
            lambda[pridx]-lambda[prrange[0]], SN[prrange[0]], unit, round(100*SN[prrange[0]]/max(SN))
      print,''

      printf, funit,'<p align="center"><table width="95%" border="0" cellspacing="0" cellpadding="0">'
      printf,funit,'<tr><td width="37%" align="left">S/N ratio at the CWL</td>'
      printf,funit,'<td width="5%" align="center">:</td><td width="58%" align="left">'
      printf,funit,format='(E10.2," ", A," (",I3," % of maximum)")', SN[cwlidx], unit, round(100*SN[cwlidx]/max(SN))
      printf,funit,'</td></tr>'
      printf,funit,'<tr><td width="37%" align="left">CWL S/N ratio at +'
      printf,funit,format='(F4.2,"A</td>")', lambda[cwlrange[1]]-lambda[cwlidx]
      printf,funit,'<td width="5%" align="center">:</td><td width="58%" align="left">'
      printf,funit,format='(E10.2," ", A," (",I3," % of maximum)")', SN[cwlrange[1]], unit, round(100*SN[cwlrange[1]]/max(SN))
      printf,funit,'</td></tr>'
      printf,funit,'<tr><td width="37%" align="left">CWL S/N ratio at -'
      printf,funit,format='(F4.2,"A</td>")', lambda[cwlidx]-lambda[cwlrange[0]]
      printf,funit,'<td width="5%" align="center">:</td><td width="58%" align="left">'
      printf,funit,format='(E10.2," ", A," (",I3," % of maximum)")', SN[cwlrange[0]], unit, round(100*SN[cwlrange[0]]/max(SN))
      printf,funit,'</td></tr>'
      printf,funit,'</table></p>'

      printf, funit,'<p align="center"><table width="95%" border="0" cellspacing="0" cellpadding="0">'
      printf,funit,'<tr><td width="37%" align="left">S/N ratio at max. pi-blue</td>'
      printf,funit,'<td width="5%" align="center">:</td><td width="58%" align="left">'
      printf,funit,format='(E10.2," ", A," (",I3," % of maximum)")', SN[pbidx], unit, round(100*SN[pbidx]/max(SN))
      printf,funit,'</td></tr>'
      printf,funit,'<tr><td width="37%" align="left">pi-blue S/N ratio at +'
      printf,funit,format='(F4.2,"A</td>")', lambda[pbrange[1]]-lambda[pbidx]
      printf,funit,'<td width="5%" align="center">:</td><td width="58%" align="left">'
      printf,funit,format='(E10.2," ", A," (",I3," % of maximum)")', SN[pbrange[1]], unit, round(100*SN[pbrange[1]]/max(SN))
      printf,funit,'</td></tr>'
      printf,funit,'<tr><td width="37%" align="left">pi-blue S/N ratio at -'
      printf,funit,format='(F4.2,"A</td>")', lambda[pbidx]-lambda[pbrange[0]]
      printf,funit,'<td width="5%" align="center">:</td><td width="58%" align="left">'
      printf,funit,format='(E10.2," ", A," (",I3," % of maximum)")', SN[pbrange[0]], unit, round(100*SN[pbrange[0]]/max(SN))
      printf,funit,'</td></tr>'
      printf,funit,'</table></p>'

      printf, funit,'<p align="center"><table width="95%" border="0" cellspacing="0" cellpadding="0">'
      printf,funit,'<tr><td width="37%" align="left">S/N ratio at max. sigma</td>'
      printf,funit,'<td width="5%" align="center">:</td><td width="58%" align="left">'
      printf,funit,format='(E10.2," ", A," (",I3," % of maximum)")', SN[sidx], unit, round(100*SN[sidx]/max(SN))
      printf,funit,'</td></tr>'
      printf,funit,'<tr><td width="37%" align="left">sigma S/N ratio at +'
      printf,funit,format='(F4.2,"A</td>")', lambda[srange[1]]-lambda[sidx]
      printf,funit,'<td width="5%" align="center">:</td><td width="58%" align="left">'
      printf,funit,format='(E10.2," ", A," (",I3," % of maximum)")', SN[srange[1]], unit, round(100*SN[srange[1]]/max(SN))
      printf,funit,'</td></tr>'
      printf,funit,'<tr><td width="37%" align="left">sigma S/N ratio at -'
      printf,funit,format='(F4.2,"A</td>")', lambda[sidx]-lambda[srange[0]]
      printf,funit,'<td width="5%" align="center">:</td><td width="58%" align="left">'
      printf,funit,format='(E10.2," ", A," (",I3," % of maximum)")', SN[srange[0]], unit, round(100*SN[srange[0]]/max(SN))
      printf,funit,'</td></tr>'
      printf,funit,'</table></p>'

      printf, funit,'<p align="center"><table width="95%" border="0" cellspacing="0" cellpadding="0">'
      printf,funit,'<tr><td width="37%" align="left">S/N ratio at max. pi-red</td>'
      printf,funit,'<td width="5%" align="center">:</td><td width="58%" align="left">'
      printf,funit,format='(E10.2," ", A," (",I3," % of maximum)")', SN[pridx], unit, round(100*SN[pridx]/max(SN))
      printf,funit,'</td></tr>'
      printf,funit,'<tr><td width="37%" align="left">pi-red S/N ratio at +'
      printf,funit,format='(F4.2,"A</td>")', lambda[prrange[1]]-lambda[pridx]
      printf,funit,'<td width="5%" align="center">:</td><td width="58%" align="left">'
      printf,funit,format='(E10.2," ", A," (",I3," % of maximum)")', SN[prrange[1]], unit, round(100*SN[prrange[1]]/max(SN))
      printf,funit,'</td></tr>'
      printf,funit,'<tr><td width="37%" align="left">pi-red S/N ratio at -'
      printf,funit,format='(F4.2,"A</td>")', lambda[pridx]-lambda[prrange[0]]
      printf,funit,'<td width="5%" align="center">:</td><td width="58%" align="left">'
      printf,funit,format='(E10.2," ", A," (",I3," % of maximum)")', SN[prrange[0]], unit, round(100*SN[prrange[0]]/max(SN))
      printf,funit,'</td></tr>'
      printf,funit,'</table></p>'
    endif

    ;----------------------------------------------------------------------------------
    ; close the device if needed
    ;----------------------------------------------------------------------------------
    if (input.spec.multiplot) then begin
      ; close the figure
      tmp = reset_graphic(figID, eps=epsfig)
      !p.multi = 0
    endif

  endfor

endif



;**********************************************************************************
;* PLOT THE PROFILES                                                              *
;**********************************************************************************

; ----------------------------------------------------------------------------------
; Loop for Sigma, Pi-red, Pi-blue, CWL pictures
; ----------------------------------------------------------------------------------
for l=0,3 do begin

  if l eq 0 then begin
    emtype   = 'CWL'         ; change names in filename
    emttl    = 'CWL'         ; change names in plot title
    emstokes = cwlstokes
  endif
  if l eq 1 then begin 
    emtype   = 'PIB'         ; change names in filename
    emttl    = 'max. pi-blue' ; change names in plot title
    emstokes = pibstokes
  endif
  if l eq 2 then begin
    emtype   = 'SIG'         ; change names in filename
    emttl    = 'max. sigma'  ; change names in plot title
    emstokes = sigstokes
  endif
  if l eq 3 then begin 
    emtype   = 'PIR'         ; change names in filename
    emttl    = 'max. pi-red' ; change names in plot title
    emstokes = pirstokes
  endif

  ;----------------------------------------------------------------------------------
  ; Only plot spectral data if requested
  ;----------------------------------------------------------------------------------
  if input.profile.intensity || input.profile.polangle || input.profile.polfrac then begin


    ;----------------------------------------------------------------------------------
    ; If multiplot is selected and some spectral plotting is requested,
    ; then we open a device for plotting here.
    ;----------------------------------------------------------------------------------
    nplots=1.
    if (input.profile.multiplot) then begin
      nplots=0.
      if input.profile.intensity  then nplots++
      if input.profile.polangle   then nplots++
      if input.profile.polfrac    then nplots++
      if input.profile.wavelength then nplots++

      ; make the figures title
      ttl = ' - Profile at '+emttl
      ttl = prefix+ttl

      ; set aspect ratio, figure zoom, fontsize and line thickness
      aspect   = 2.5/sqrt(nplots*3.)
      figzoom  = 1.2
      fontsize = 12*sqrt(nplots)
      lth      = 2.0
	 
      ; open a device for plotting the profiles as function of R
      if (input.general.plot_print eq 1) then begin	
        figID = init_graphic(aspect=aspect, fontsize=fontsize, fontname=fontname, zoom=figzoom, title=ttl)
      endif
      if (input.general.plot_print eq 2) then begin	
        fname = '_profiles'+emtype+'.eps'
        fname = pictdir+'/'+psprefix+fname
        figID = init_graphic(aspect=aspect, fontsize=fontsize, fontname=fontname, /PSfont, zoom=figzoom, eps=fname, bits_per_pixel=8)
      endif
      !p.multi=[0,1,nplots]
      !p.background=truecolor('white')	; use a white background
      !p.color     =truecolor('black')	; use a black foreground
      !x.margin    = [10,9]
    endif

    ;----------------------------------------------------------------------------------
    ; Convert the Stokes vectors in intensity, polarised fraction,
    ; polarisation angle and norm. S/N 
    ;----------------------------------------------------------------------------------
    intem   = reform(emstokes[0,*])
    pfracem = reform(sqrt(emstokes[1,*]^2+emstokes[2,*]^2)/emstokes[0,*])
    SNem    = reform(sqrt((emstokes[1,*]^2+emstokes[2,*]^2)/emstokes[0,*]))
    gammaem = 0.5*reform(atan(emstokes[2,*],emstokes[1,*]))*!radeg
    EM      = emstokes[4,*]

    ;----------------------------------------------------------------------------------
    ; print some information to the screen and the html-file
    ;----------------------------------------------------------------------------------
    print,''
    print,format='("  - Profile information at ",A)',emttl
    printf, funit,format='("<h2>Profile information at ",A,"</h2>")',emttl
    printf, funit,'<p align="center"><table width="95%" border="0" cellspacing="0" cellpadding="0">'
    printf,funit,'<tr><td><b>Channel ID</b></td><td><b>R (m)</b></td><td><b>Wavelength (A)</b></td>'
    if filterflag eq 0 then begin
      print,'    Channel ID  |  R (m)  | Wavelength (A) | Intensity (ph/s/A) | Pol. angle (deg.) | Pol. fraction'
      printf,funit,'<td><b>Intensity (ph/s/A)</b></td><td><b>Pol. angle (deg.)</b></td><td><b>Pol. fraction</b></td></tr>'
    endif else begin
      print,'    Channel ID  |  R (m)  | Wavelength (A) |  Intensity (ph/s)  | Pol. angle (deg.) | Pol. fraction'
      printf,funit,'<td><b>Intensity (ph/s)</b></td><td><b>Pol. angle (deg.)</b></td><td><b>Pol. fraction</b></td></tr>'
    endelse
    print,'    -----------------------------------------------------------------------------------------'
    for k =0,nchan-1 do begin
      print,format='("      ",A7, "   |  ",F6.3," |     ",F7.2,"    |     ",E10.1,"     |     ",F7.2,"       |  ",F5.1,"%")',$
            chanID[k], R_res[2,k], EM[k], intem[k],gammaem[k], pfracem[k]*100.
      printf,funit,format='("<tr><td>",A,"</td><td>",F6.3,"</td><td>",F7.2,"</b></td>")',$
            chanID[k], R_res[2,k], EM[k]
       printf,funit,format='("<td>",E10.1,"</td><td>",F7.2,"</td><td>",F5.1,"%</td></tr>")',$
            intem[k],gammaem[k], pfracem[k]*100.
    endfor
    printf,funit,'</table></p>'


    ;----------------------------------------------------------------------------------
    ; plot intensity profile if requested 
    ;----------------------------------------------------------------------------------
    if input.profile.intensity then begin

      ;----------------------------------------------------------------------------------
      ; make the figures title
      ;----------------------------------------------------------------------------------
      ttl = ' - Intensity at '+emttl
      ttl = prefix+ttl

      ;----------------------------------------------------------------------------------
      ; If multiplot is NOT selected, then we open a device for plotting here.
      ;----------------------------------------------------------------------------------
      if ~(input.profile.multiplot) then begin
        ; set aspect ratio, figure zoom, fontsize and line thickness
        aspect   = 1.5
        figzoom  = 1.2
        fontsize = 12
        lth      = 2.0
		 
        ; open a device for plotting the profile as function of R
        if (input.general.plot_print eq 1) then begin	
          figID = init_graphic(aspect=aspect, fontsize=fontsize, fontname=fontname, zoom=figzoom)
        endif
        if (input.general.plot_print eq 2) then begin	
          fname = '_intensity'+emtype+'.eps'
          fname = pictdir+'/'+psprefix+fname
          figID = init_graphic(aspect=aspect, fontsize=fontsize, fontname=fontname, /PSfont, zoom=figzoom, eps=fname, bits_per_pixel=8)
        endif
        !p.multi     = 0
        !p.background= truecolor('white')	; use a white background
        !p.color     = truecolor('black')	; use a black foreground
      endif

      ;----------------------------------------------------------------------------------
      ; the plotting
      ;----------------------------------------------------------------------------------

      ; set limits
      xmin = min(R)
      xmax = max(R)
      ymin = 0.
      ymax = 1.05*max(intem)
      ; plot intensity profile
      if filterflag eq 0 then yttl='Intensity [photons/s/Angstrom]' else yttl='Intensity [photons/s]'
      plot, R_res[2,*], intem, thick=lth, psym=-1, $
            xs=1, xr=[xmin,xmax], ys=1, yr=[ymin,ymax],$
            xtitle='R - Centre of Mass [m]', ytitle=yttl, title=ttl

      ;----------------------------------------------------------------------------------
      ; close the device if needed
      ;----------------------------------------------------------------------------------
      if ~(input.profile.multiplot) then begin
        ; close the figure
        tmp = reset_graphic(figID, eps=epsfig)
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
      ttl = ' - Pol. Angle at '+emttl+' / Residual angle'
      ttl = prefix+ttl

      ;----------------------------------------------------------------------------------
      ; If multiplot is NOT selected, then we open a device for plotting here.
      ;----------------------------------------------------------------------------------
      if ~(input.profile.multiplot) then begin
        ; set aspect ratio, figure zoom, fontsize and line thickness
        aspect   = 1.5
        figzoom  = 1.2
        fontsize = 12
        lth      = 2.0
		 
        ; open a device for plotting the profile as function of R
        if (input.general.plot_print eq 1) then begin	
          figID = init_graphic(aspect=aspect, fontsize=fontsize, fontname=fontname, zoom=figzoom)
        endif
        if (input.general.plot_print eq 2) then begin	
          fname = '_polangle'+emtype+'.eps'
          fname = pictdir+'/'+psprefix+fname
          figID = init_graphic(aspect=aspect, fontsize=fontsize, fontname=fontname, /PSfont, zoom=figzoom, eps=fname, bits_per_pixel=8)
        endif
        !p.multi     = 0
        !p.background= truecolor('white')	; use a white background
        !p.color     = truecolor('black')	; use a black foreground
      endif


      ;----------------------------------------------------------------------------------
      ; the plotting of the pol. angle profile
      ;----------------------------------------------------------------------------------
      ; set limits
      xmin = min(R)
      xmax = max(R)
      ymin = min(gammaem)-0.05*(max(gammaem)-min(gammaem))
      ymax = max(gammaem)+0.05*(max(gammaem)-min(gammaem))
      if (ymax-ymin) lt 2. then begin
        ymax=mean(gammaem)+1.
        ymin=mean(gammaem)-1.
      endif
      ; plot pol. angle profile
      plot, R_res[2,*], gammaem, thick=lth, psym=-1,color=truecolor('black'),$
            xs=1, xr=[xmin,xmax], ys=5, yr=[ymin,ymax],$
            xtitle='R - Centre of Mass [m]',$
            title=ttl
      axis, yaxis=0, ys=1, yr=[ymin,ymax], ytitle='Pol. angle '+textoidl('\gamma')+' [degrees]'
      ; oplot ideal polarization profile 0, -90, -180 (because alpha0 defined between 0 and !pi)
      oplot, R_res[2,*], [alpha0/!dtor,alpha0/!dtor],         color=truecolor('black'), thick=0.5, linestyle=3
      oplot, R_res[2,*], [alpha0/!dtor-90,alpha0/!dtor-90],   color=truecolor('black'), thick=0.5, linestyle=3
      oplot, R_res[2,*], [alpha0/!dtor-180,alpha0/!dtor-180], color=truecolor('black'), thick=0.5, linestyle=3
      ; Residual
      position = [!x.window, !y.window]
      position = position[[0,2,1,3]]
      resid    = gammaem-alpha0*!radeg
      ymin     = -1
      ymax     =  1
      plot, R_res[2,*], resid, linestyle=2, thick=lth, psym=-1,     $
            xs=5, xr=[xmin,xmax], ys=5, yr=[ymin,ymax],color=truecolor('red'),     $
            position=position, /noerase
      oplot,R_res[2,*], resid+90, color=truecolor('red'), thick=lth, linestyle=2, psym=-1
      oplot,R_res[2,*], resid+180, color=truecolor('red'), thick=lth, linestyle=2, psym=-1
      oplot, !x.crange, [0,0], col=truecolor('black'), linestyle=1, thick=0.5*lth
      axis, yaxis=1, ys=1, yr=[ymin,ymax], ytitle='Residual angle [degrees]'
      ; legend
      makelegend, ['full-simulation polarisation angle','ideal polarisation angle','residual angle'],                                        $
               textcolors=truecolor(['black','black','red']),                                       $
               colors    =truecolor(['black','black','red']),                                       $
               linestyle=[0,3,2],thick=[lth,0.5,lth],                                             $
               charsize=1./sqrt(nplots), spacing=1./sqrt(nplots), /clear, /box, /top, /left
   
      ;----------------------------------------------------------------------------------
      ; close the device if needed
      ;----------------------------------------------------------------------------------
      if ~(input.profile.multiplot) then begin
        ; close the figure
        tmp = reset_graphic(figID, eps=epsfig)
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
      ttl = ' - Pol. fraction/SN-ratio at '+emttl
      ttl = prefix+ttl

      ;----------------------------------------------------------------------------------
      ; If multiplot is NOT selected, then we open a device for plotting here.
      ;----------------------------------------------------------------------------------
      if ~(input.profile.multiplot) then begin
        ; set aspect ratio, figure zoom, fontsize and line thickness
        aspect   = 1.5
        figzoom  = 1.2
        fontsize = 12
        lth      = 2.0
		 
        ; open a device for plotting the profile as function of R
        if (input.general.plot_print eq 1) then begin	
          figID = init_graphic(aspect=aspect, fontsize=fontsize, fontname=fontname, zoom=figzoom)
        endif
        if (input.general.plot_print eq 2) then begin	
          fname = '_polfrac'+emtype+'.eps'
          fname = pictdir+'/'+psprefix+fname
          figID = init_graphic(aspect=aspect, fontsize=fontsize, fontname=fontname, /PSfont, zoom=figzoom, eps=fname, bits_per_pixel=8)
        endif
        !p.multi     = 0
        !p.background= truecolor('white')	; use a white background
        !p.color     = truecolor('black')	; use a black foreground
        !x.margin    = [10,9]
      endif


      ;----------------------------------------------------------------------------------
      ; the plotting of the pol. fraction (only where the polarised intensity is nonzero)
      ;----------------------------------------------------------------------------------
      ; set limits
      xmin = min(R)
      xmax = max(R)
      ymin = 0.0
      ymax = 1.05
      ; polarised fraction
      plot, R_res[2,*], pfracem, thick=lth, psym=-1,$
            xs=1, xr=[xmin,xmax], ys=5, yr=[ymin,ymax],$
            xtitle='R - Centre of Mass [m]', title=ttl
      oplot, !x.crange, [1,1], col=truecolor('black'), linestyle=1, thick=0.5*lth
      axis, yaxis=0, ys=1, yr=[ymin,ymax], ytitle='Linear polarisation fraction'
      ; S/N ratio
      position = [!x.window, !y.window]
      position = position[[0,2,1,3]]
      ymin=0
      ymax=1.05*max(SNem)
      plot, R_res[2,*], SNem, linestyle=2, thick=lth, psym=-1, $
            xs=5, xr=[xmin,xmax], ys=5, yr=[ymin,ymax],     $
            position=position, /noerase
      axis, yaxis=1, ys=1, yr=[ymin,ymax], ytitle='S/N-ratio [(photons/s/Angstrom)^{1/2}]'
      ; legend
      makelegend, ['linear polarisation fraction', 'S/N ratio'],    $
              textcolors=truecolor(['black','black']),          $
              colors    =truecolor(['black','black']),          $
              linestyle=[0,2],thick=[lth,lth],                  $
              charsize=1./sqrt(nplots), spacing=1./sqrt(nplots),$
              /clear, /box, /top, /left

      ;----------------------------------------------------------------------------------
      ; close the device if needed
      ;----------------------------------------------------------------------------------
      if ~(input.profile.multiplot) then begin
        ; close the figure
        tmp = reset_graphic(figID, eps=epsfig)
        !p.multi = 0
      endif

    endif

    ;----------------------------------------------------------------------------------
    ; plot wavelength profile if requested 
    ;----------------------------------------------------------------------------------
    if input.profile.wavelength then begin

      ;----------------------------------------------------------------------------------
      ; make the figures title
      ;----------------------------------------------------------------------------------
      ttl = ' - Wavelength at '+emttl
      ttl = prefix+ttl

      ;----------------------------------------------------------------------------------
      ; If multiplot is NOT selected, then we open a device for plotting here.
      ;----------------------------------------------------------------------------------
      if ~(input.profile.multiplot) then begin
        ; set aspect ratio, figure zoom, fontsize and line thickness
        aspect   = 1.5
        figzoom  = 1.2
        fontsize = 12
        lth      = 2.0
		 
        ; open a device for plotting the profile as function of R
        if (input.general.plot_print eq 1) then begin	
          figID = init_graphic(aspect=aspect, fontsize=fontsize, fontname=fontname, zoom=figzoom)
        endif
        if (input.general.plot_print eq 2) then begin	
          fname = '_wavelen'+emtype+'.eps'
          fname = pictdir+'/'+psprefix+fname
          figID = init_graphic(aspect=aspect, fontsize=fontsize, fontname=fontname, /PSfont, zoom=figzoom, eps=fname, bits_per_pixel=8)
        endif
        !p.multi     = 0
        !p.background= truecolor('white')	; use a white background
        !p.color     = truecolor('black')	; use a black foreground
        !x.margin    = [10,9]
      endif


      ;----------------------------------------------------------------------------------
      ; the plotting of the pol. fraction (only where the polarised intensity is nonzero)
      ;----------------------------------------------------------------------------------
      ; set limits
      xmin = min(R)
      xmax = max(R)
      ymin = min(EM) -0.05*(max(EM)-min(EM))
      ymax = max(EM) +0.05*(max(EM)-min(EM))
      if (ymax-ymin) lt 2. then begin
        ymin = mean(EM)-1.
        ymax = mean(EM)+1.
      endif
      ; plot wavelength profile
      if filterflag eq 0 then yttl='Wavelength [Angstrom]' else yttl='Filter central wavelength [Angstrom]'
      plot, R_res[2,*], EM, thick=lth, psym=-1,$
            xs=1, xr=[xmin,xmax], ys=1, yr=[ymin,ymax],$
            xtitle='R - Centre of Mass [m]', title=ttl, ytitle=yttl

      ;----------------------------------------------------------------------------------
      ; close the device if needed
      ;----------------------------------------------------------------------------------
      if ~(input.profile.multiplot) then begin
        ; close the figure
        tmp = reset_graphic(figID, eps=epsfig)
        !p.multi = 0
      endif

    endif

    ;----------------------------------------------------------------------------------
    ; close the device if needed
    ;----------------------------------------------------------------------------------
    if (input.profile.multiplot) then begin
      ; close the figure
      tmp = reset_graphic(figID, eps=epsfig)
      !p.multi = 0
    endif

  endif
endfor

;----------------------------------------------------------------------------------
; close the html text-file that stored the screen output
;----------------------------------------------------------------------------------
; close the txtfile:
printf, funit, '</body>'
printf, funit, '</html>'
free_lun, funit
close, funit

;----------------------------------------------------------------------------------
; Reset to default plotting
;----------------------------------------------------------------------------------
!p.multi = 0

end
