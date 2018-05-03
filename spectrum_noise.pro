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

pro spectrum_noise, settingfile, inputfile, outputfile

; routine calculates the noise and the associated statistical errors in the 
; polarisation angles
;
; v1.0, mdebock 23/08/2007
;

;----------------------------------------------------------------------------------
; Distinguish between Windows and UNIX directory separators
;----------------------------------------------------------------------------------
if strcmp( !version.os,'Win32',/fold_case) then sep='\' else sep='/'

;**********************************************************************************
;* READ FILTER SETTINGS AND INPUTFILE                                             *
;**********************************************************************************

;----------------------------------------------------------------------------------
; Physical constants
;----------------------------------------------------------------------------------
input    = read_setting('..'+sep+'physics'+sep+'physic_constants.xml')
c        = input.physics.c		; light speed (m/s)
e        = input.physics.e		; electron charge (C)
h        = input.physics.h		; Planck constant (Js)

;----------------------------------------------------------------------------------
; Transmision and detector settings
;----------------------------------------------------------------------------------
input = read_setting(settingfile)

ts          = input.measurementsettings.sampletime
pem_op      = input.measurementsettings.PEMoperationpoint

t_PEM       = input.transmission.PEM
t_polariser = input.transmission.polariser
t_lens      = input.transmission.lens
t_fibre     = input.transmission.fibre
t_filter    = input.transmission.filterscope

d_QE        = input.detector.QuantumEfficiency
d_exnoise   = input.detector.ExcessNoise	; decreases the effective quantum efficiency
d_cap       = input.detector.capacity
d_amp       = input.detector.amplification
d_darkcurr  = input.detector.darkcurrent


;----------------------------------------------------------------------------------
; Read the inputdata
;----------------------------------------------------------------------------------
restore, inputfile

; get some inportant data
nchan   = n_elements(tspec[0,*])				; number of channels
dlambda = lambda[1]-lambda[0]					; step in wavelength
nlambda = n_elements(lambda)					; number of wavelengths

; we need the find the filter intensity of the total spectrum, the filtered
; polarised fraction, the polarisation angle and the wavelength at maximum polarised 
; intensity for pi-blue, sigma and pi-red. And this for every channel.
Int        = fltarr(nchan,3)
Frac       = fltarr(nchan,3)	; the '3' stands for the blue shifted pi, the sigma and the red shifted pi
Angle      = fltarr(nchan,3)
Wavelength = fltarr(nchan,3)
for k=0,nchan-1 do begin
  ;----------------------------------------------------------------------------------
  ; Find out the 'nonzero' polarised intensities
  ;----------------------------------------------------------------------------------
  tpnonzero = where(tpspec[*,k] gt 0.0)

  ;----------------------------------------------------------------------------------
  ; Find the polarisation fraction and intesities
  ;----------------------------------------------------------------------------------
  ; polarisation fraction 
  tpfrac = fltarr(nlambda)
  tpfrac[tpnonzero] = tpspec[tpnonzero,k]/tspec[tpnonzero,k]
  pmins  = ppspec[*,k] - spspec[*,k]

  ; polarised intensity = (pol. fraction)*sqrt(total intenstity)
  tpint = fltarr(nlambda)
  tpint[tpnonzero] = tpspec[tpnonzero,k]/sqrt(tspec[tpnonzero,k])

  ;----------------------------------------------------------------------------------
  ; Find the positions of max. pi-blue, max. sigma and max.pi-red
  ;----------------------------------------------------------------------------------
  ; first guess where the maximum intensities are
  mx = min(pmins,sidx)		; approximate position of the sigma peak
  mx = max(pmins[0:sidx], pbidx)	; approximate position of the blue pi peak
  mx = max(pmins[sidx:*], pridx)	; approximate position of the red pi peak
  pridx = pridx+sidx
  range = round(0.5/dlambda)	; we take a range of +- 0.5 Angstrom (i.e. 1 Angstrom in total)

  ; then find actual positions of the maxima
  startidx = pbidx-range
  stopidx  = pbidx+range
  mx = max(tpint[startidx:stopidx], pbidx)
  pbidx = startidx + pbidx

  startidx = sidx-range
  stopidx  = sidx+range
  mx = max(tpint[startidx:stopidx], sidx)
  sidx = startidx + sidx

  startidx = pridx-range
  stopidx  = pridx+range
  mx = max(tpint[startidx:stopidx], pridx)
  pridx = startidx + pridx

  ;----------------------------------------------------------------------------------
  ; Finally save the total intensity, pol. fraction, pol. angle
  ; and the wavelength at these positions
  ;----------------------------------------------------------------------------------
  Int[k,0]   = tspec[pbidx,k]
  Frac[k,0]  = tpfrac[pbidx]
  Angle[k,0] = alpha[pbidx,k]

  Int[k,1]    = tspec[sidx,k]
  Frac[k,1]   = tpfrac[sidx]
  Angle[k,1]  = alpha[sidx,k]

  Int[k,2]    = tspec[pridx,k]
  Frac[k,2]   = tpfrac[pridx]
  Angle[k,2]  = alpha[pridx,k]

  Wavelength[k,0] = lambda[pbidx]
  Wavelength[k,1] = lambda[sidx]
  Wavelength[k,2] = lambda[pridx]

endfor

;-----------------------------------------------------------------------------------------------------
; Now find out what there is to know about the photons that fall onto the detector
;-----------------------------------------------------------------------------------------------------
photons_sec = t_PEM * t_polariser * t_lens * t_fibre * t_filter * Int	; in photons/second
intensity   = photons_sec * h*c/(Wavelength*1e-10)			; in Watt
photon_curr = e * photons_sec * d_QE 			; the current in the detector due to the photons
photon_volt = d_amp* photon_curr * ts/d_cap				; the voltage signal of the detector is the current
									; integrated over the sampletime and divided by the capacity
									; and multiplied by the amplification

;-----------------------------------------------------------------------------------------------------
; Find out the noise
;-----------------------------------------------------------------------------------------------------
; the noise of the detector is defined by equation (1-5) in the document:
;   SD-28 - Characteristics and use of Si APD (Avalanche Photodiode) - by Hamamatsu
; it is:
;  noise_curr^2 = 2*e*(photon_curr + dark_curr[0])*d_amp^2*d_exnoise*bandwidth + 2*e*dark_curr[1]*bandwidth
;
;

; apart from the noise due to the photon statistics, a part of the noise comes from the detector dark current.
; this dark current can be seen as caused by 'dark' photons.
dark_photons = (ts/e) * d_darkcurr
; the 'number' of (Poisson) noise is 'photons' then is
noise_photons = sqrt(photons + dark_photons)
; the noise voltage is then given by the noise charge divided by the capacity and multiplied by the amplification
noise_volt = d_amp * e * (d_QE/d_exnoise) * noise_photons / d_cap

;-----------------------------------------------------------------------------------------------------
; The signal and the noise 
;-----------------------------------------------------------------------------------------------------
; The signals of interest to us (i.e. the ones needed to get the polarisation angle) are the frequency
; components 2*omega, where omega is the modulation frequency PEM1 and PEM2, respectively.
; These signals are proportional to the polarised fraction of the emission, the second Besel function
; of the PEM's operation point, and to the cosine and sine of twice the polarisation angle respectively:
I2w1 = photon_volt * Frac * beselj(pem_op,2)/sqrt(2.0) * abs(cos(2*Angle))
I2w2 = photon_volt * Frac * beselj(pem_op,2)/sqrt(2.0) * abs(sin(2*Angle))

; this means the 2 signal to noise ratios do depend on the polarisation angle.
; However, based upon the derivation found the JET Report 'JET-R(96)10' by Nick Hawkes we find
; that the noise/error in the polarisation angle does not depend on the polarisation angle:
alpha_noise = 1.0 / (beselj(pem_op,2) * sqrt(2.0)) * noise_volt/(Frac*photon_volt)
; (see 'Stark Simulation Code' report for the full derivation of this statistical error)


;-----------------------------------------------------------------------------------------------------
; Lets plot 
;-----------------------------------------------------------------------------------------------------
pseudocol
!p.background=9

window, 0, xsize=1000, ysize=800

xmin = min(R_res[2,*]) - 0.05*(max(R_res[2,*])-min(R_res[2,*]))
xmax = max(R_res[2,*]) + 0.05*(max(R_res[2,*])-min(R_res[2,*]))
ymin = 0.0
ymax = 1e8

plot, [xmin,xmax],[ymin,ymax], color=0, xs=1, ys=1, /nodata,$
      xtitle='R (m)', ytitle='# photons', charsize=2.5, charthick=2
oplot, R_res[2,*], photons[*,0], color=3, thick=2, linestyle=0
oplot, R_res[2,*], photons[*,1], color=0, thick=2, linestyle=0
oplot, R_res[2,*], photons[*,2], color=1, thick=2, linestyle=0
;legend
xyouts, xmin+0.05*(xmax-xmin), ymin+0.9*(ymax-ymin),charsize=2.5,charthick=2,color=3, textoidl('\pi') + '- blue'
xyouts, xmin+0.05*(xmax-xmin), ymin+0.85*(ymax-ymin),charsize=2.5,charthick=2,color=0, textoidl('\sigma')
xyouts, xmin+0.05*(xmax-xmin), ymin+0.8*(ymax-ymin),charsize=2.5,charthick=2,color=1, textoidl('\pi') + '- red'


window, 1, xsize=1000, ysize=800
ymin = 0.0
ymax = 1e3

plot, [xmin,xmax],[ymin,ymax], color=0, xs=1, ys=1, /nodata,$
      xtitle='R (m)', ytitle='signal/noise', charsize=2.5, charthick=2
oplot, R_res[2,*], I2w1[*,0]/noise_volt[*,0], color=3, thick=2, linestyle=0
oplot, R_res[2,*], I2w2[*,0]/noise_volt[*,0], color=3, thick=2, linestyle=2
oplot, R_res[2,*], I2w1[*,1]/noise_volt[*,0], color=0, thick=2, linestyle=0
oplot, R_res[2,*], I2w2[*,1]/noise_volt[*,0], color=0, thick=2, linestyle=2
oplot, R_res[2,*], I2w1[*,2]/noise_volt[*,0], color=1, thick=2, linestyle=0
oplot, R_res[2,*], I2w2[*,2]/noise_volt[*,0], color=1, thick=2, linestyle=2
;legend
xyouts, xmin+0.05*(xmax-xmin), ymin+0.9*(ymax-ymin),charsize=2.5,charthick=2,color=3, textoidl('\pi, I_{2 \omega 1}') + '-  blue'
xyouts, xmin+0.05*(xmax-xmin), ymin+0.85*(ymax-ymin),charsize=2.5,charthick=2,color=3, textoidl('\pi, I_{2 \omega 2}') + '-- blue'
xyouts, xmin+0.05*(xmax-xmin), ymin+0.8*(ymax-ymin),charsize=2.5,charthick=2,color=0, textoidl('-  \sigma, I_{2 \omega 1}')
xyouts, xmin+0.05*(xmax-xmin), ymin+0.75*(ymax-ymin),charsize=2.5,charthick=2,color=0, textoidl('-- \sigma, I_{2 \omega 2}')
xyouts, xmin+0.05*(xmax-xmin), ymin+0.7*(ymax-ymin),charsize=2.5,charthick=2,color=1, textoidl('-  red \pi, I_{2 \omega 1}')
xyouts, xmin+0.05*(xmax-xmin), ymin+0.65*(ymax-ymin),charsize=2.5,charthick=2,color=1, textoidl('-- red \pi, I_{2 \omega 2}')


window, 2, xsize=1000, ysize=800
ymin = 0.0
ymax = 1.0

plot, [xmin,xmax],[ymin,ymax], color=0, xs=1, ys=1, /nodata,$
      xtitle='R (m)', ytitle='Error in pol. angle' + textoidl('\alpha'), charsize=2.5, charthick=2
oplot, R_res[2,*], alpha_noise[*,0]*!radeg, color=3, thick=2, linestyle=0
oplot, R_res[2,*], alpha_noise[*,1]*!radeg, color=0, thick=2, linestyle=0
oplot, R_res[2,*], alpha_noise[*,2]*!radeg, color=1, thick=2, linestyle=0
;legend
xyouts, xmin+0.05*(xmax-xmin), ymin+0.9*(ymax-ymin),charsize=2.5,charthick=2,color=3, textoidl('\pi') + '- blue'
xyouts, xmin+0.05*(xmax-xmin), ymin+0.85*(ymax-ymin),charsize=2.5,charthick=2,color=0, textoidl('- \sigma')
xyouts, xmin+0.05*(xmax-xmin), ymin+0.8*(ymax-ymin),charsize=2.5,charthick=2,color=1, textoidl('\pi') + 'red'



;/////////////////////////////
; HIC SUNT LEONES
;/////////////////////////////




;-----------------------------------------------------------------------------------------------------
; Save the result
;-----------------------------------------------------------------------------------------------------
;save, xyz0, B_v0, B_xyz, B_w, B_vec, C_xyz, C_vec, vec0,$	; central points data
;      Bfld0, Efld0, Dshift0, Sshift0, alpha0, psi0,$
;      gp_xyz, gp_vel, gp_vec, gp_Bfld, gp_Efld,$		; grid point data
;      gp_alpha, gp_psi, gp_emis,$
;      R, Z, psi, RZpsi, Rm, RZ_emis, R_emis, psi_emis,$		; resolution data
;      psi_res, R_res,$
;      lambda, tspec, tpspec, ppspec, puspec,$			; spectral data
;      spspec, suspec, alpha, palpha, salpha,$
;      filename=outputfile


end
