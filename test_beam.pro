pro test_beam, calc=calc
; tokamak and field parameters
R0    = 0.83
a     = 0.60
Bphi  = -0.5
q0    = 1.0
qa    = 3.0
qidx  = 3.0
Bp0   = 0.0
Bpa   = 0.0
Bpidx = 1.0
shafr = 0.1
elong = 1.66

; beam parameters
B0   =[0,-2,0]
w0   =0.1
xi   =70
delta=90
chi  =1.2

Bdens0=5e15
edens0=5e19
Qion  =5e-20
Qemit =2.5e-12

Brange =[0.4,1.90]
nl     = 20
nw     = 20
ntheta = 18

dl     = (Brange[1]-Brange[0])/(nl-1)
dw     = (3*w0)/(nw-1)
dtheta = 2*!pi/ntheta

;beamfile = 'beam/9168swbes60kevpini.xdr'
beamfile = 'beam/beam_bes_PINI_S_18501_t0.29s_63keV.xdr'
readbeam = not(keyword_set(calc))

; degrees to radians
d2r = !pi/180

; beam vector
Bv = [cos(d2r*xi)*sin(d2r*delta), sin(d2r*xi)*sin(d2r*delta), cos(d2r*delta)]
X2 = coordtrans([1,0,0],[[0,0,0],[Bv]])


; calculate or readthe beam
npt  = long(nl)*long(nw)*long(ntheta)
Bpt  =fltarr(3,npt)
m=0
for i=0,nl-1 do begin
for j=0,nw-1 do begin
for k=0,ntheta-1 do begin
   theta = k*dtheta 
   w   = j*dw
   l   = Brange[0] + i*dl

   Bwv = coordtrans([0,w*cos(theta),w*sin(theta)],[[0,0,0],[X2]])
   Bpt[*,m] = B0 + l*Bv + Bwv

   m++
endfor
endfor
endfor

if readbeam then begin
   beam = read_beam(Bpt,Bv,beamfile)
endif else begin
   beam = calc_beam(Bpt, B0,Bv,d2r*chi,w0,$
                    Bdens0,edens0,Qion,Qemit,$
                    R0, a, shafr, elong,$
                    Bphi,q0,qa,qidx,Bp0,Bpa,Bpidx)
endelse
; Beam density, electron density and emission rate as function of Bpt:
Bdens = beam.Bdens
edens = beam.edens
Qemit = beam.emission
; The same as function of length along, distance from and angle around the beam axis
Bdens3 =fltarr(nl,nw,ntheta)
edens3 =fltarr(nl,nw,ntheta)
Qemit3 =fltarr(nl,nw,ntheta)
m=0
for i=0,nl-1 do begin
for j=0,nw-1 do begin
for k=0,ntheta-1 do begin
   Bdens3[i,j,k] = Bdens[m]
   edens3[i,j,k] = edens[m]
   Qemit3[i,j,k] = Qemit[m]
   m++
endfor
endfor
endfor


; Beam density plot
;------------------
; plot the geometry
  pseudocol,/on
  !p.background=9
  !p.multi=[0,2,2]
  window, 0, xsize=1000, ysize=1000,$					; open a window for plotting the MSE geometry
          title='Beam density'
  plot_3Dbox, [0.0,2.0], [0.0,0.0], [0.0,0.0],$			; plot the X-axis and titles
              color=0,thick=3,$
              az=-60,ax=40,$
              xs=1,xr=[-0.1,2.1],$
              ys=1,yr=[-2.1,0.51],$
              zs=1,zr=[-0.6,0.6], $
              charsize=2.0, charthick=1.2,$
              xtitle='X (m)', ytitle='Y (m)', ztitle='Z (m)'
  plots, [0.0,0.0], [-2.0,0.0], [0.0,0.0], color=0, thick=3, /T3D	; plot the Y-axis
  plots, [0.0,0.0], [0.0,0.0], [-0.5,0.5], color=0, thick=3, /T3D	; plot the Z-axis

; plot the beam density in 3D
  Bdensmax = max(Bdens)
  Bdensmin = min(Bdens)
  if abs(Bdensmax-Bdensmin) gt 1e-6 then cfactor  = 245/(Bdensmax-Bdensmin)$
  else cfactor = 245/(2*Bdensmax)
  for m=0,npt-1 do begin
    c = round(10+ cfactor*(Bdens[m]-Bdensmin))
    plots, Bpt[0,m],Bpt[1,m],Bpt[2,m],$
           color=c, psym=4, symsize=2, /T3D
  endfor


  ; plot the beam density in the equatorial plane (Z=0)
  l  = Brange[0]+findgen(nl)*dl
  d  = findgen(nw)*dw
  d2 = findgen(2*nw-1)*dw - (nw-1)*dw
  lddens =fltarr(nl,2*nw-1)
  idxpi =ntheta/2
  lddens[0:nl-1,0:nw-1]    = rotate(Bdens3[*,*,idxpi[0]],7)
  lddens[0:nl-1,nw:2*nw-2] = Bdens3[*,1:*,0]
  surface, lddens ,l, d2, color=0,charsize=2.0, charthick=1.2,$
           az=-30,xtitle='Length along the beam (m)',$
           zs=1,zr=[0,1.1*Bdensmax],$
           ytitle='Distance from beam axis (m)', ztitle='Beam density (m^-3)'

  ; plot the central beam density as a funtion of lenght along the beam
   ldens = Bdens3[*,0,0]
   plot, l, ldens, color=0, thick=3, charsize=1.6,$
        xs=1, xr=[0,1.1*Brange[1]],ys=1,yr=[0,1.1*Bdensmax],$
        xtitle='Length along the beam (m)', ytitle='Beam density (m^-3)'

  ; plot the beam density profile at (max.) 6 position along the beam
   if nl lt 9 then step=1 else step=ceil(nl/6) 
   ls = 0
  ; plot the axes
    plot,[0,2*w0],[0,1.1*Bdensmax], /nodata, color=0, charsize=1.6,$
        xs=1, xr=[0,3*w0],ys=1,yr=[0,1.1*Bdensmax],$
        xtitle='Distance from beam axis (m)', ytitle='Beam density (m^-3)'
  ; plot the data
   for i=0,nl-1,step do begin
     ddens = Bdens3[i,*,0]
     oplot, d, ddens, color=0, linestyle=ls, thick=3
     ls=ls+1
   endfor
  ; plot the legend
   str = string(FORMAT='("-   L = ",(F4.2)," m")',l[0])
   xyouts, 2.0*w0, 0.9*Bdensmax,str,color=0
   if nl gt 1 then begin
     str = string(FORMAT='("..   L = ",(F4.2)," m")',l[1*step])
     xyouts, 2.0*w0, 0.8*Bdensmax,str,color=0
   endif
   if nl gt 2 then begin
     str = string(FORMAT='("--  L = ",(F4.2)," m")',l[2*step])
     xyouts, 2.0*w0, 0.7*Bdensmax,str,color=0
   endif
   if nl gt 3 then begin
     str = string(FORMAT='("-.   L = ",(F4.2)," m")',l[3*step])
     xyouts, 2.0*w0, 0.6*Bdensmax,str,color=0
   endif
   if nl gt 4 then begin
     str = string(FORMAT='("-..  L = ",(F4.2)," m")',l[4*step])
     xyouts, 2.0*w0, 0.5*Bdensmax,str,color=0
   endif
   if nl gt 4 then begin
     str = string(FORMAT='("_ _ L = ",(F4.2)," m")',l[5*step])
     xyouts, 2.0*w0, 0.4*Bdensmax,str,color=0
   endif

; electron density plot
;----------------------
; plot the geometry
  !p.background=9
  !p.multi=[0,2,2]
  window, 1, xsize=1000, ysize=1000,$					; open a window for plotting the MSE geometry
          title='electron density'
  plot_3Dbox, [0.0,2.0], [0.0,0.0], [0.0,0.0],$			; plot the X-axis and titles
              color=0,thick=3,$
              az=-60,ax=40,$
              xs=1,xr=[-0.1,2.1],$
              ys=1,yr=[-2.1,0.1],$
              zs=1,zr=[-0.6,0.6], $
              charsize=2.0, charthick=1.2,$
              xtitle='X (m)', ytitle='Y (m)', ztitle='Z (m)'
  plots, [0.0,0.0], [-2.0,0.0], [0.0,0.0], color=0, thick=3, /T3D	; plot the Y-axis
  plots, [0.0,0.0], [0.0,0.0], [-0.5,0.5], color=0, thick=3, /T3D	; plot the Z-axis

; plot the electron density in 3D
  edensmax = max(edens)
  edensmin = min(edens)
  if abs(edensmax-edensmin) gt 1e-6 then cfactor  = 245/(edensmax-edensmin)$
  else cfactor = 245/(2*edensmax)
  for m=0,npt-1 do begin
    c = round(10+ cfactor*(edens[m]-edensmin))
    plots, Bpt[0,m],Bpt[1,m],Bpt[2,m],$
           color=c, psym=4, symsize=2, /T3D
  endfor

  ; plot the electron density in the equatorial plane (Z=0)
  l  = Brange[0]+findgen(nl)*dl
  d  = findgen(nw)*dw
  d2 = findgen(2*nw-1)*dw - (nw-1)*dw
  lddens =fltarr(nl,2*nw-1)
  idxpi =ntheta/2
  lddens[0:nl-1,0:nw-1]    = rotate(edens3[*,*,idxpi[0]],7)
  lddens[0:nl-1,nw:2*nw-2] = edens3[*,1:*,0]
  surface, lddens ,l, d2, color=0,charsize=2.0, charthick=1.2,$
           az=-30,xtitle='Length along the beam (m)',$
           zs=1,zr=[0,1.1*edensmax],$
           ytitle='Distance from beam axis (m)', ztitle='electron density (m^-3)'

  ; plot the central electron density as a funtion of lenght along the beam
   ldens = edens3[*,0,0]
   plot, l, ldens, color=0, thick=3, charsize=1.6,$
        xs=1, xr=[0,1.1*Brange[1]],ys=1,yr=[0,1.1*edensmax],$
        xtitle='Length along the beam (m)', ytitle='electron density (m^-3)'

  ; plot the electron density profile at (max.) 6 position along the beam
   if nl lt 9 then step=1 else step=ceil(nl/6) 
   ls = 0
  ; plot the axes
    plot,[0,2*w0],[0,1.1*edensmax], /nodata, color=0, charsize=1.6,$
        xs=1, xr=[0,3*w0],ys=1,yr=[0,1.1*edensmax],$
        xtitle='Distance from beam axis (m)', ytitle='electron density (m^-3)'
  ; plot the data
   for i=0,nl-1,step do begin
     ddens = edens3[i,*,0]
     oplot, d, ddens, color=0, linestyle=ls, thick=3
     ls=ls+1
   endfor
  ; plot the legend
   str = string(FORMAT='("-   L = ",(F4.2)," m")',l[0])
   xyouts, 2.0*w0, 0.9*edensmax,str,color=0
   if nl gt 1 then begin
     str = string(FORMAT='("..   L = ",(F4.2)," m")',l[1*step])
     xyouts, 2.0*w0, 0.8*edensmax,str,color=0
   endif
   if nl gt 2 then begin
     str = string(FORMAT='("--  L = ",(F4.2)," m")',l[2*step])
     xyouts, 2.0*w0, 0.7*edensmax,str,color=0
   endif
   if nl gt 3 then begin
     str = string(FORMAT='("-.   L = ",(F4.2)," m")',l[3*step])
     xyouts, 2.0*w0, 0.6*edensmax,str,color=0
   endif
   if nl gt 4 then begin
     str = string(FORMAT='("-..  L = ",(F4.2)," m")',l[4*step])
     xyouts, 2.0*w0, 0.5*edensmax,str,color=0
   endif
   if nl gt 4 then begin
     str = string(FORMAT='("_ _ L = ",(F4.2)," m")',l[5*step])
     xyouts, 2.0*w0, 0.4*edensmax,str,color=0
   endif


; Emission rate plot
;-------------------
; plot the geometry
  !p.background=9
  !p.multi=[0,2,2]
  window, 2, xsize=1000, ysize=1000,$					; open a window for plotting the MSE geometry
          title='Emission rate'
  plot_3Dbox, [0.0,2.0], [0.0,0.0], [0.0,0.0],$			; plot the X-axis and titles
              color=0,thick=3,$
              az=-60,ax=40,$
              xs=1,xr=[-0.1,2.1],$
              ys=1,yr=[-2.1,0.1],$
              zs=1,zr=[-0.6,0.6], $
              charsize=2.0, charthick=1.2,$
              xtitle='X (m)', ytitle='Y (m)', ztitle='Z (m)'
  plots, [0.0,0.0], [-2.0,0.0], [0.0,0.0], color=0, thick=3, /T3D	; plot the Y-axis
  plots, [0.0,0.0], [0.0,0.0], [-0.5,0.5], color=0, thick=3, /T3D	; plot the Z-axis

; plot the emission rate
  Qemitmax = max(Qemit)
  Qemitmin = min(Qemit)
  if abs(Qemitmax-Qemitmin) gt 1e-6 then cfactor  = 245/(Qemitmax-Qemitmin)$
  else cfactor = 245/(2*Qemitmax)
  for m=0,npt-1 do begin
    c = round(10+ cfactor*(Qemit[m]-Qemitmin))
    plots, Bpt[0,m],Bpt[1,m],Bpt[2,m],$
           color=c, psym=4, symsize=2, /T3D
  endfor

  ; plot the emission rate in the equatorial plane (Z=0)
  l  = Brange[0]+findgen(nl)*dl
  d  = findgen(nw)*dw
  d2 = findgen(2*nw-1)*dw - (nw-1)*dw
  ldemit =fltarr(nl,2*nw-1)
  idxpi =ntheta/2
  ldemit[0:nl-1,0:nw-1]    = rotate(Qemit3[*,*,idxpi[0]],7)
  ldemit[0:nl-1,nw:2*nw-2] = Qemit3[*,1:*,0]
  surface, ldemit ,l, d2, color=0,charsize=2.0, charthick=1.2,$
           az=-30,xtitle='Length along the beam (m)',$
           zs=1,zr=[0,1.1*Qemitmax],$
           ytitle='Distance from beam axis (m)', ztitle='Emission rate (photons/m^3/s)'

  ; plot the central emission rate as a funtion of lenght along the beam
   lemit = Qemit3[*,0,0]
   plot, l, lemit, color=0, thick=3, charsize=1.6,$
        xs=1, xr=[0,1.1*Brange[1]],ys=1,yr=[0,1.1*Qemitmax],$
        xtitle='Length along the beam (m)', ytitle='Emission rate (photons/m^3/s)'

  ; plot the emission rate profile at (max.) 6 position along the beam
   if nl lt 9 then step=1 else step=ceil(nl/6)
   ls = 0
  ; plot the axes
    plot,[0,2*w0],[0,1.1*Qemitmax], /nodata, color=0, charsize=1.6,$
        xs=1, xr=[0,3*w0],ys=1,yr=[0,1.1*Qemitmax],$
        xtitle='Distance from beam axis (m)', ytitle='Emission rate (photons/m^3/s)'
  ; plot the data
   for i=0,nl-1,step do begin
     demit = Qemit3[i,*,0]
     oplot, d, demit, color=0, linestyle=ls, thick=3
     ls=ls+1
   endfor
  ; plot the legend
   str = string(FORMAT='("-   L = ",(F4.2)," m")',l[0])
   xyouts, 2.0*w0, 0.9*Qemitmax,str,color=0
   if nl gt 1 then begin
     str = string(FORMAT='("..   L = ",(F4.2)," m")',l[1*step])
     xyouts, 2.0*w0, 0.8*Qemitmax,str,color=0
   endif
   if nl gt 2 then begin
     str = string(FORMAT='("--  L = ",(F4.2)," m")',l[2*step])
     xyouts, 2.0*w0, 0.7*Qemitmax,str,color=0
   endif
   if nl gt 3 then begin
     str = string(FORMAT='("-.   L = ",(F4.2)," m")',l[3*step])
     xyouts, 2.0*w0, 0.6*Qemitmax,str,color=0
   endif
   if nl gt 4 then begin
     str = string(FORMAT='("-..  L = ",(F4.2)," m")',l[4*step])
     xyouts, 2.0*w0, 0.5*Qemitmax,str,color=0
   endif
   if nl gt 4 then begin
     str = string(FORMAT='("_ _ L = ",(F4.2)," m")',l[5*step])
     xyouts, 2.0*w0, 0.4*Qemitmax,str,color=0
   endif

!p.multi=0
end