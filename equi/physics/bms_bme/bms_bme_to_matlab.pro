pro bms_bme_to_matlab
;+
; BME_BMS_TO_MATLAB
;
; Routine exports beam stopping and beam emission rates (<sigma.v> in m^3/s) to a MATLAB m-file
;-

  ; Load the IDL saved file
  restore, file='bms_bme.sav', /v
  ; contains:
  ;  density    : [n_density x 1]     range of densities    [m^{-3}]
  ;  temperature: [n_temperature x 1] range of temperatures [eV]
  ;  energy     : [n_energy x 1]      range of energies/amu [eV]
  ;  sigmav_bms : [n_energy x n_density x n_temperature] beam stopping rate
  ;               as function of density temperature and energy [m^3/s]
  ;  sigmav_bme : [n_energy x n_density x n_temperature] beam emission rate
  ;               as function of density temperature and energy [m^3/s]
 
  ; write to matlab
  wr2matlab, density, 'density', 'bms_bme'
  wr2matlab, temperature, 'temperature', 'bms_bme', /append
  wr2matlab, energy, 'energy', 'bms_bme', /append
  wr2matlab, sigmav_bms, 'sigmav_bms', 'bms_bme', /append
  wr2matlab, sigmav_bme, 'sigmav_bme', 'bms_bme', /append

end
