pro get_bms_bme_adas, Erange, Nrange, Trange, zeff=zeff, fraction=fraction, filename=filename
;+
; GET_BMS_BME_ADAS, Erange, Trange, Nrange, zeff=zeff, fraction=fraction, filename=filename
; 
; Gets the beam stopping (bms) and beam emission (bes) rate coefficients (<sigma v>) from the
; ADAS data base for a defined ranges of beam energies, plasma densities and temperatures.
; Zeff or the fractions of H/D, C and O are optional keywords.
; The result is saved as an IDL file with name 'filename'. This file will contain the variables
;   energy     : type=double,  1D array
;               The beam energy/amu vector [eV]
;   density    : type=double,  1D array
;               The plasma density vector [m^-3]
;   temperature: type=double,  1D array
;               The plasma temperatures [eV]
;   sigmav_bms : type=double,  3D array
;               [n_energy x n_density x n_temperature] array with <sigma v> for the beam stopping [m^3/s]
;   sigmav_bme : type=double,  3D array
;               [n_energy x n_density x n_temperature] array with <sigma v> for the beam emission [m^3/s]
;
; :Input:
;   Erange     : required, type=double, 1D array
;               [Emin, Emax, #elements, lin/log] 4 -element array determining the energy vector [eV].
;               e.g. [1e4, 1e5, 20, 0] will generate a 20 element energy vector, linearly increasing
;               form 1e4eV to 1e5eV
;   Nrange     : required, type=double, 1D array
;               [Nmin, Nmax, #elements, lin/log] 4 -element array determining the density vector [m^{-3}].
;               e.g. [1e18, 1e20, 19, 1] will generate a 20 element density vector, exponentially increasing
;               form 1e18m^{-3} to 1e20m^{-3}
;   Trange     : required, type=double, 1D array
;               [Tmin, Tmax, #elements, lin/log] 4 -element array determining the temperature vector [eV].
;               e.g. [10, 10000, 19, 1] will generate a 20 element temperature vector, exponentially increasing
;               form 10eV to 1e4eV
; :Keywords:
;   zeff       : optional, type=double, scalar/1D array
;               Zeff, based on carbon and oxygen.
;               If zeff is a scalar, then n_O = 0.5*n_C is assumed (cf.: R. Neu, K. Asmussen et al.
;                "Monitor for the carbon and oxygen impurities in the ASDEX Upgrade tokamak"
;                Rev. Sci. Instrum. 67, 1829 (1996); doi:10.1063/1.1147522)
;               If zeff is a 2-element array than the first element is Zeff and the second element
;               is the ratio n_O/n_C
;               If zeff is omitted then a Zeff=1 is assumed
;  fraction    : optional, type=double, 1D array
;               The fractions of H/D, C and O. e.g. [0.94, 0.06, 0.02] means n_D = 0.94*n_i, n_C = 0.06*n_i, n_O = 0.02*n_i
;               and n_e = n_D + 6*n_C + 8*n_O = n_i * (0.94 + 0.06*6 + 0.02*8)
;               The fraction keyword overrules the zeff keyword. If omitted, fraction is calculated from zeff
;  filename    : required, type=string, scalar
;               The name of the file that the results will be saved to.
; :History:
;   03/05/2011 - v1.0 - Maarten De Bock (m.f.m.d.bock@tue.nl)
;                Creation of the routine. WARNING: not working! => find out how to access ADAS first
;-

  fraction = [1.0, 0.0, 0.0]
  restore, f='bms_h.dat', /verbose      

  ; Beam emission
  filesbme  = [ '/funsrv1/home/adas/adas/adf22/bme97#h/bme97#h_h1.dat', $
                '/funsrv1/home/adas/adas/adf22/bme97#h/bme97#h_c6.dat', $
                '/funsrv1/home/adas/adas/adf22/bme97#h/bme97#h_o8.dat' ]
  ; Beam stopping
  filesbms  = [ '/funsrv1/home/adas/adas/adf21/bms97#h/bms97#h_h1.dat', $
                '/funsrv1/home/adas/adas/adf21/bms97#h/bms97#h_c6.dat', $
               '/funsrv1/home/adas/adas/adf21/bms97#h/bms97#h_o8.dat' ]

  SIGMAV=fltarr(n_elements(energy),n_elements(density) , n_elements(temperature))       
  SIGMAVbes=fltarr(n_elements(energy),n_elements(density) , n_elements(temperature))  
  for zz=0,n_elements(density)-1 do begin
      dens=fltarr(n_elements(energy))+density[zz]
       for ff=0,n_elements(temperature)-1 do  begin  
         te=fltarr(n_elements(energy))+temperature[ff]
                    
;               
               read_adf21, files=filesbms, energy=energy, te=te,  dens=dens,   $
               fraction=fraction, data=data 
               sigmav[*,zz,ff]=data

                read_adf22, files=filesbme, energy=energy, te=te,  dens=dens,   $
                fraction=fraction,data=data  
                sigmavbes[*,zz,ff]=data

          endfor
    print,zz
    endfor

save, f='bms_test.dat',energy,density,temperature,sigmav,sigmavbes,/xdr
stop


end




