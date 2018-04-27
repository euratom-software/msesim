;*****************************************************
;The version for unix with some comments
;uses ADAS data for BS cross-sections
;uses int_tabulated function instead of qromb
;*****************************************************
pro grid,shot, save=save, infile=infile,noplasmacalc=noplasmacalc,bes=bes, mastu=mastu

;-------------------------------------------------------------------------------
; some initialisation for the plotting
;-------------------------------------------------------------------------------
pseudocol	; use my personal pseudocolor table:
!p.background=9 ;  * 0...9   : black,red,green,blue,cyaan,yellow,magenta,dark gray,light gray,white
                ;  * 10...253: shades from white, over yellow, red and magenta to blue
                ;  * 254,255 : white

;-------------------------------------------------------------------------------
; Initialise the common bloc that will contains the variables used by all
; subroutines
;-------------------------------------------------------------------------------
@com_bloc

;-------------------------------------------------------------------------------
; Some settings (that are not in the input file)
;-------------------------------------------------------------------------------
Enumber=3	; number of energy components of the beam (this should always be 3 I guess: full, half, third)
plotE=0		; which energy component should be plotted (0=full, 1=half, 2=third)
tcalc=[.29]     ; Times at which we want to calculate the beam
ORNL_PINI=2	; Oakridge Beam (1) or PINI-beam (2)
cnt_surf=1	; if 1 => contour plots are made, if 2 => surface plots are made 

; Following 4 settings are only used in a 'no plasma' calculation.
; (in case of plasma they are read from the shotdata)
abeam=2		  ; Beam mass number (H=1, D=2)
beamvoltage=6.3e4 ; Beam voltage (V)
beampower=[1500.,150.,45.0]*1e3   ; Beam power in full, half and third energy components (W)
beamtime=[0.0]	  ; time vector of the beam voltage
                  ; (arbitrary for the 'no plasma' calculation, just set it to 0.0)
; again this is just for the 'no plasma' calculation (in case of plasma this is overwritten)
currentfull =beampower[0]/beamvoltage
currenthalf =beampower[1]/(beamvoltage/2.0)
currentthird=beampower[2]/(beamvoltage/3.0)

;-------------------------------------------------------------------------------
; Read the settings file
;-------------------------------------------------------------------------------
print, ''
print,'Reading setting file'
print,'--------------------'
readdat, infile=infile

;-------------------------------------------------------------------------------
; initialise the output arrays and the xc,yc,zc arrays
;-------------------------------------------------------------------------------
noplasma = fltarr(xdata[2], ydata[2], zdata[2])
i_space = fltarr(xdata[2], ydata[2], zdata[2], n_elements(tcalc),Enumber)

step_x = abs(xdata[1]-xdata[0])/((xdata[2]-1)>1)
step_y = abs(ydata[1]-ydata[0])/((ydata[2]-1)>1)
step_z = abs(zdata[1]-zdata[0])/((zdata[2]-1)>1)
xc = xdata[0]+findgen(xdata[2])*step_x
yc = ydata[0]+findgen(ydata[2])*step_y
zc = zdata[0]+findgen(zdata[2])*step_z

;-------------------------------------------------------------------------------
; Set the type of beam
;-------------------------------------------------------------------------------
if ORNL_PINI eq 1 then setbeam $	; set Oakridge beam
else begin
  if keyword_set(mastu) then setpini_mastu, X0=mastu else setpini      ; set PINI beam
endelse

;-------------------------------------------------------------------------------
; Read the shot data (only if 'no plasma' is NOT set)
;-------------------------------------------------------------------------------
if ~keyword_set(noplasmacalc) then begin
  print, ''
  print,'Reading shot data'
  print,'-----------------'
  ; read in the shot data (Te, ne, beam mass, beam voltage and beam currents)
  readshot,shot
  print, ''
  print,'Rescaling cross sections'
  print,'------------------------'
  print,''
  ; scale the cross sections
  rescale_cross_sec
endif

;-------------------------------------------------------------------------------
; The main loop of calculation
;-------------------------------------------------------------------------------
print, ''
print,'Calculating Beam'
print,'----------------'

windowcnt=4	; for plotting: there are already 4 plots made (1 in setpini, 3 in readshot)
; loop through the calculation times
for tind=0,n_elements(tcalc)-1 do begin
  ; we do not use the actual calculation time, but
  ; the time in the beam voltage time vector that comes
  ; closest to the calculation time
  tmp=min(abs(beamtime-tcalc[tind]),time)
  print,format='("* time        : ",F5.3," s")',beamtime[time]

  ; The (usually 3) beam energy components at this time are:
  ebeam=beamvoltage[time]/(findgen(Enumber)+1)
  print,format='("  Beam energy : ",F4.1," keV (full),  ",F4.1," keV (half),  ",F4.1," keV (third)")',$
        ebeam[0]*1e-3, ebeam[1]*1e-3, ebeam[2]*1e-3

  ; The full beam energy should be at least 5 keV
  if (ebeam[0] le 5e3) then begin
    i_space(*, *, *, tind, *) = 0.
    print,'  WARNING: Full beam energy should be higher then 5 keV. Skipping calculation for this time point!'
    ; skip to the next iteration of the for-loop
    continue	
  end

  ; loop through the xc, yc and zc grid points
  print, format='($,"  Calculating beam at each x,y,z-gridpoint ")'
  for i=0,xdata[2]-1 do begin
    x = xc[i]
    for j=0,ydata(2)-1 do begin
      y = yc[j]
      for k=0,zdata(2)-1 do begin
        z = zc[k]

        ; calculate the solid angles of the beam particles
        solidangle,x,y,z,solid,solidsingle
        ; the spread of the beam is given by this solid angle:
        noplasma[i, j, k]=solid[0]

        if keyword_set(noplasmacalc) then begin 
          ; in case of no plasma, the beam current distribution is only determined by this spread
          i_space[i, j, k, tind, *] = solid[0]
        endif else begin 
          ; in case of plasma, the beam current distribution is determined by this spread
          ; and the fact that beam particles get ionised when they move through the plasma.
          ; so we need to calculate the beam attenuation:
          integral, x, y, z, solidsingle, int
          i_space(i, j, k, tind, *) = int
        endelse
      endfor
    endfor
    print, format='($,".")'
  endfor
  print, ''


  ; do some plotting
  ;-----------------
  ; some contour and color settings
  ncnt   = 100					; number of contours
  colors = round(findgen(ncnt)/(ncnt-1)*243)+10	; colors from white over yellow, red and magenta to blue

  ; make contour plots at 5 z-positions
  ;------------------------------------
  ; distribute the 5 z-position equally over the zc-range
  nz   = n_elements(zc)
  zind = round(findgen(5)/4*(nz-1))
  ; open a new window
  window,windowcnt, xsize=1300, ysize=500,$
         title=string(format='("Beam current density at 5 z-positions for t=",F5.3,"s")',beamtime[time])
  windowcnt++
  !p.multi=[0,5,1]
  ; we'll use the same maximum for all 5 plots (so you can see the attenuation more easily)
  maxdata = max(i_space[*,*,zind,tind,plotE])
  for i=0,4 do begin
    contourdata=i_space[*,*,zind[i],tind,plotE]
    if cnt_surf then begin
      contour,contourdata,xc,yc,/iso,/fill,nlevel=ncnt,c_colors=colors,charsize=2,$
              zr=[0, maxdata],zs=1,color=0,xtitle='x (cm)', ytitle='y (cm)',$
              title='z='+string(format='(I3)',zc[zind[i]])+' cm'
    endif else begin
      surface,contourdata,xc,yc,charsize=4.0,color=0,$
              zr=[0, maxdata],zs=1,xtitle='x (cm)', ytitle='y (cm)',$
              title='z='+string(format='(I3)',zc[zind[i]])+' cm'
    endelse
  endfor

  ; make contours at plot at 6 y-positions
  ;---------------------------------------
  ; distribute the 6 y-position equally over the yc-range
  ny   = n_elements(yc)
  yind = round(findgen(6)/5*(ny-1))
  ; open a new window
  window,windowcnt, xsize=1300, ysize=700,$
         title=string(format='("Beam current density at 6 y-positions for t=",F5.3,"s")',beamtime[time])
  windowcnt++
  !p.multi=[0,3,2]
  ; we'll use the same maximum for all 6 plots (so you can see the attenuation more easily)
  maxdata = max(i_space[*,yind,*,tind,plotE])
  for i=0,5 do begin
    contourdata=transpose(i_space[*,yind[i],*,tind,plotE],[0,2,1])
    if cnt_surf then begin
      contour,contourdata,xc,zc,/iso,/fill,nlevel=ncnt,c_colors=colors,charsize=2,$
              zr=[0, maxdata],color=0,zs=1,xtitle='x (cm)', ytitle='z (cm)',$
              title='y='+string(format='(I4)',yc[yind[i]])+' cm'
    endif else begin
      surface,contourdata,xc,zc,charsize=6.0,color=0,$
              xs=1,ys=1,zr=[0, maxdata],zs=1,xtitle='x (cm)', ytitle='z (cm)',$
              title='y='+string(format='(I4)',yc[yind[i]])+' cm'
    endelse
  endfor
endfor


;-------------------------------------------------------------------------------
; convert beam current density into beam neutral density
;-------------------------------------------------------------------------------
print, ''
print,'Converting beam current density into beam neutral density'
print,'---------------------------------------------------------'

; do the conversion
neutr_den_beam,i_space,neutral
; in case 'no plasma' was NOT selected: map Te and ne onto the xc,yc,zc grid
if ~keyword_set(noplasmacalc) then rescale_nt_new,nel_space,tel_space
print,'done!'

;-------------------------------------------------------------------------------
; calculate the fraction of the beam density that is in the n=2 excited state
;-------------------------------------------------------------------------------
if ~keyword_set(noplasmacalc) then begin   ; if there is no plasma then there is no
  print, ''                                ; excitation into the n=2 state by electron collisions!
  print,'Calculate n=2-population '
  print,'------------------------'
  n2,nel_space,tel_space,n2
  print,'done!'
endif

;-------------------------------------------------------------------------------
; If the keyword 'bes' is set: calculate the emission rate
;-------------------------------------------------------------------------------
if keyword_set(bes) then begin
  print, ''
  print,'Calculating the beam emission rate'
  print,'----------------------------------'

  ; calculating the emission rate is only possible when there is plasma:
  if keyword_set(noplasmacalc) then begin
     print,'  WARNING: The "noplasmacalc" keyword is set!'
     print,'           Therefore the beam emission rate can not be calculated!'
  endif else begin
    ; calculate beam emission rate in photons/s/m^3
    bes,nel_space,tel_space,ratebes
    ratebes=ratebes*neutral/1e6
    print, 'done!'
    ; do some plotting
    ;-----------------
    ; some contour and color settings
    ncnt   = 100					; number of contours
    colors = round(findgen(ncnt)/(ncnt-1)*243)+10	; colors from white over yellow, red and magenta to blue

    ; make a plot for every calculated time
    ;--------------------------------------
    for t=0,n_elements(tcalc)-1 do begin
      ; we do not use the actual calculation time, but
      ; the time in the beam voltage time vector that comes
      ; closest to the calculation time
      tmp=min(abs(beamtime-tcalc[t]),time)

      ; make contour plots at 5 z-positions
      ;------------------------------------
      ; distribute the 5 z-position equally over the zc-range
      nz   = n_elements(zc)
      zind = round(findgen(5)/4*(nz-1))
      ; open a new window
      window,windowcnt, xsize=1300, ysize=500,$
            title=string(format='("Beam emission rate at 5 z-positions for t=",F5.3,"s")',beamtime[time])
      windowcnt++
      !p.multi=[0,5,1]
      ; we'll use the same maximum for all 5 plots (so you can see the attenuation more easily)
      maxdata = max(ratebes[*,*,zind,t,plotE])
      for i=0,4 do begin
        contourdata=ratebes[*,*,zind[i],t,plotE]
        if cnt_surf then begin
          contour,contourdata,xc,yc,/iso,/fill,nlevel=ncnt,c_colors=colors,charsize=2,$
                  zr=[0, maxdata],zs=1,color=0,xtitle='x (cm)', ytitle='y (cm)',$
                  title='z='+string(format='(I3)',zc[zind[i]])+' cm'
        endif else begin
          surface,contourdata,xc,yc,charsize=4.0,color=0,$
                  zr=[0, maxdata],zs=1,xtitle='x (cm)', ytitle='y (cm)',$
                  ztitle='Beam emission rate (photons/s/m^3)',$
                  title='z='+string(format='(I3)',zc[zind[i]])+' cm'
        endelse
      endfor

      ; make contours at plot at 6 y-positions
      ;---------------------------------------
      ; distribute the 6 y-position equally over the yc-range
      ny   = n_elements(yc)
      yind = round(findgen(6)/5*(ny-1))
      ; open a new window
      window,windowcnt, xsize=1300, ysize=700,$
            title=string(format='("Beam emission rate at 6 y-positions for t=",F5.3,"s")',beamtime[time])
      windowcnt++
      !p.multi=[0,3,2]
      ; we'll use the same maximum for all 6 plots (so you can see the attenuation more easily)
      maxdata = max(ratebes[*,yind,*,t,plotE])
      for i=0,5 do begin
        contourdata=transpose(ratebes[*,yind[i],*,t,plotE],[0,2,1])
        if cnt_surf then begin
          contour,contourdata,xc,zc,/iso,/fill,nlevel=ncnt,c_colors=colors,charsize=2,$
                  zr=[0, maxdata],color=0,zs=1,xtitle='x (cm)', ytitle='z (cm)',$
                  title='y='+string(format='(I4)',yc[yind[i]])+' cm'
        endif else begin
          surface,contourdata,xc,zc,charsize=6.0,color=0,$
                  xs=1,ys=1,zr=[0, maxdata],zs=1,xtitle='x (cm)', ytitle='z (cm)',$
                  ztitle='Beam emission rate (photons/s/m^3)',$
                  title='y='+string(format='(I4)',yc[yind[i]])+' cm'
        endelse
      endfor
    endfor
  endelse
endif

;-------------------------------------------------------------------------------
; If the keyword 'bes' is set: calculate the emission rate
;-------------------------------------------------------------------------------
if keyword_set(save) then begin
  print, ''
  print,'Save to xdr-file'
  print,'----------------'

  ; create a default name based upon the type of beam and calculation
  if keyword_set(noplasmacalc) then begin
    shotstr = 'noplasma'
  endif else begin
    shotstr     = strcompress(string(shot),/remove_all)
  endelse
  if ORNL_PINI then begin
    if whatbeam eq 'ss' then beamtype='ORNL_S'
    if whatbeam eq 'sw' then beamtype='ORNL_SW'
  endif else begin
    if whatbeam eq 'ss' then beamtype='PINI_S'
    if whatbeam eq 'sw' then beamtype='PINI_SW'
  endelse
  if keyword_set(bes) && ~keyword_set(noplasmacalc) then begin
    calctype='bes'
  endif else begin
    calctype='dens'
  endelse
  defaultname = 'beam_'+calctype+'_'+beamtype+'_'+shotstr+'.xdr'

  ; offer the user the possibility to change this default name:
  fname = ''
  read,prompt='Give a filename (default="'+defaultname+'") : ',fname
  fname = strtrim(fname,2)
  if strcmp(fname,'') then fname = defaultname

  ; save the file
  if keyword_set(bes) && ~keyword_set(noplasmacalc) then begin
    save, i_space,shot,whatbeam,neutral,n2,nel_space,tel_space,xc,yc,zc,$
          beamtime,nel,tel,rshot,zshot,currentfull,currenthalf,currentthird,$
          beamvoltage,abeam, ratebes,$
          filename=fname
  endif else begin
    if keyword_set(noplasmacalc) then begin
      save, i_space,shot,whatbeam,neutral,xc,yc,zc,$
            beamtime,currentfull,currenthalf,currentthird,$
            beamvoltage,abeam,$
            filename=fname
    endif else begin
      save, i_space,shot,whatbeam,neutral,n2,nel_space,tel_space,xc,yc,zc,$
            beamtime,nel,tel,rshot,zshot,currentfull,currenthalf,currentthird,$
            beamvoltage,abeam,$
            filename=fname
    endelse
  endelse

endif
print,''
;-------------------------------------------------------------------------------
; reset plotting parameters to default ones
;-------------------------------------------------------------------------------
pseudocol,/off	; switch off my personal pseudocolor table
!p.background=0 ; set the background to black
!p.multi=0	; disable multiple plots in one window

end




